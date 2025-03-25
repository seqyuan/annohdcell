# annohdcell

Python package for processing HD spatial transcriptomics data

## Installation

```bash
pip install annohdcell
```

## Usage

### Command Line Interface

Process HD spatial data from bins to cells:
```bash
annohdcell bin2cell \
  -d /path/to/square_002um \
  -t /path/to/image.tif \
  -o /path/to/output
```

Transform h5ad file for R compatibility:
```bash
annohdcell h5trans \
  -c /path/to/input.h5ad \
  -o /path/to/output.h5ad
```

## Documentation

Full API documentation is available at [ReadTheDocs](https://annohdcell.readthedocs.io)

## License

MIT
