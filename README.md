# organage

A package to estimate organ-specific biological age using SomaScan plasma proteomics data

## System requirements

### Hardware requirements

organage package requires only a standard computer with enough RAM to support the in-memory operations.

### Software requirements

#### OS Requirements

This package is supported for macOS and Linux. The package has been tested on the following systems:
- maxOS 11.7.1
- Linux: CentOS 7.x

#### Python dependencies

- python>=3.9
- dill>=0.3.6
- pandas>=1.5.3
- scikit-learn==1.0.2. aging models were trained using this specific version of scikit-learn

## Installation

```bash
$ pip install organage
```

## Usage

- see docs/example.ipynb

## License

`organage` was created by Hamilton Oh. It is licensed under the terms of the MIT license.

## Credits

`organage` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
