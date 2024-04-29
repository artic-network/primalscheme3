from collections.abc import Iterator
import time

from tqdm import tqdm


class ProgressTracker(tqdm):
    def __init__(self, parent, process, iterable, *args, **kwargs):
        self.parent = parent
        self.process = process
        super().__init__(iterable, *args, **kwargs)

    def __iter__(self) -> Iterator:
        obj = super().__iter__()

        for i in obj:
            self.parent.signal()
            yield i


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

    def create_sub_progress(self, iter, process, *args, **kwargs) -> ProgressTracker:
        """Create a progress tracker"""
        self._subprocess = ProgressTracker(
            self, *args, iterable=iter, process=process, **kwargs, desc=process
        )
        return self._subprocess


if __name__ == "__main__":
    pm = ProgressManager()

    print(pm.process(), pm.n(), pm.total())

    for _ in pm.create_sub_progress(iter=range(10), process="test"):
        print(_)
        print(pm.process(), pm.n(), pm.total())
        time.sleep(0.1)

    print(pm.n())
