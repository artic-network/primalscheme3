from collections.abc import Iterator

from tqdm import tqdm


class ProgressTracker(tqdm):
    def __init__(self, parent, process, iterable, chrom, *args, **kwargs):
        self.parent = parent
        self.process = process
        self.chrom = chrom
        super().__init__(iterable, *args, **kwargs)

    def __iter__(self) -> Iterator:
        obj = super().__iter__()

        for i in obj:
            self.parent.signal()
            yield i

    def manual_update(
        self,
        n: int | None = None,
        total: int | None = None,
        process: str | None = None,
        chrom: str | None = None,
    ):
        # Update the progress bar manually
        if n is not None:
            self.n = n
        if total is not None:
            self.total = total
        if process is not None:
            self.process = process
        if chrom is not None:
            self.chrom = chrom

        self.parent.signal()


class ProgressManager:
    _subprocess: None | ProgressTracker

    def signal(self):
        pass

    def __init__(self):
        self._status = None
        self._subprocess = None

    def n(self) -> int | None:
        if self._subprocess:
            return self._subprocess.n
        return None

    def total(self) -> int | None:
        if self._subprocess:
            return self._subprocess.total
        return None

    def process(self) -> str | None:
        if self._subprocess:
            return self._subprocess.process
        return None

    def chrom(self) -> str | None:
        if self._subprocess:
            return self._subprocess.chrom
        return None

    def create_sub_progress(
        self, iter, chrom, process, *args, **kwargs
    ) -> ProgressTracker:
        """Create a progress tracker"""
        self._subprocess = ProgressTracker(
            self,
            *args,
            iterable=iter,
            chrom=chrom,
            process=process,
            **kwargs,
            desc=process,
        )
        return self._subprocess
