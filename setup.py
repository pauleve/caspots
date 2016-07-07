
from setuptools import setup, find_packages

version = "0.1"

setup(name="caspots",
    version = version,
    description = "Boolean network identification from time-series data",
    long_description = "",
    install_requires = [
#        "gringo",
        "caspo",
    ],
    packages = find_packages(),
    include_package_data = True,
    entry_points = {
        "console_scripts": [
            "caspots=caspots.console:run",
        ]
    }
)


