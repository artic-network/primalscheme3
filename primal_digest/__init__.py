import tomllib

# Read in some values from pyproject.toml
with open("./pyproject.toml", "rb") as toml:
    data = tomllib.load(toml)
    __version__ = data["tool"]["poetry"]["version"]
