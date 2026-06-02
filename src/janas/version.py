# File: version.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology


# src/janas/version.py
try:
    # Python >= 3.8
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    # older Python
    from pkg_resources import get_distribution, DistributionNotFound
    def version(package):
        try:
            return get_distribution(package).version
        except DistributionNotFound:
            raise PackageNotFoundError

def get_version():
    try:
        return version("janas")
    except PackageNotFoundError:
        return "0.0.0"

