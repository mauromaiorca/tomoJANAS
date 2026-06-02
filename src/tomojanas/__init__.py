# tomoJANAS — sub-volume extraction and analysis framework for Electron Tomography.
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology

try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:  # pragma: no cover - very old Python
    from pkg_resources import get_distribution, DistributionNotFound

    def version(package):  # type: ignore
        try:
            return get_distribution(package).version
        except DistributionNotFound:
            raise PackageNotFoundError


def get_version() -> str:
    try:
        return version("tomojanas")
    except PackageNotFoundError:
        return "0.0.0"


__version__ = get_version()
