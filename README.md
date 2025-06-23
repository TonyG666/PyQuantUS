# Adding a New Parser in `utc_loaders`

**Target Audience**: Software Developers with basic ultrasound knowledge

## Overview

The `utc_loaders` directory is organized so that each subfolder corresponds to a specific ultrasound probe or data format (e.g., `canon_iq`, `ge_rf`). Each folder contains the code necessary to parse and load data from that probe into the PyQuantUS analysis pipeline.

This modular approach allows the software to support multiple probe types (linear, curvilinear, etc.) and makes it easy to add support for new data formats.

---

## 🚀 Steps to Add a New Parser

### 1. Create a New Folder

- In `utc_loaders`, create a new folder named after your probe or data format (e.g., `ge_rf`).


### 2. Mimic the Structure of an Existing Parser

- Identify the model/type of the ultrasound probe (e.g., curvilinear/convex or linear).
  - _Example_: GE Logiq E10 system, using a **convex** probe
- If the probe is curvilinear/convex, **scan conversion** codes are needed.
- Use an existing folder as a template:
  - For probes requiring scan conversion, reference the structure of the `canon_iq` parser.

Typical contents of a parser folder:
```
<parser_name>/
├── __pycache__/
├── main.py          # class used to parse data from specific ultrasound machines
├── objects.py       # data structures (usually stable across parsers)
├── parser.py        # main parser module
```


### 3. Modify Contents in the Parser Folder

#### `main.py`
- Adjust for differences in the number of instance attributes depending on whether scan conversion is needed.

#### `parser.py`
- In the function that reads file information:
  - Hardcode the frequencies of the data based on internet research of the scanning device.
  - Add scan conversion codes if needed.
  - Fully understand the sample data and how to locate it programmatically.
- In the function that reads file images:
  - Add scan conversion codes if necessary.
- **Phantom Data**:
  - Even if phantom data isn't available, always write code for parsing phantom data.


### 4. Test the Parser

- Locate the sample data file for your probe.
- Open `utc_demo.ipynb` in the `CLI-Demos` folder.
- Replace:
  - `scan_path` with the path to your new data
  - `phantom_path` with the phantom data path (or reuse `scan_path` if phantom data is missing)
- Run the analysis pipeline using your parser to ensure everything works correctly.


### 5. Document the Parser

- Add or update **docstrings** for your functions and classes.
- Describe:
  - The data format
  - Any special considerations (e.g., missing metadata, unconventional layout, etc.)

---

## 📁 Example Folder Structure

```
pyquantus/
├── ...
└── image_loading/
    └── utc_loaders/
        ├── canon_iq/
        │   ├── __pycache__/
        │   ├── parser.py
        │   ├── objects.py
        │   └── main.py
        └── ge_rf/
            ├── __pycache__/
            ├── parser.py
            ├── objects.py
            └── main.py
```

---

## 🔑 Key Points

- **Consistency**: Follow the interface and conventions of existing parsers.
- **Original vs Phantom Data**: Always parse both original and phantom data for analysis.
- **Extensibility**: Design code so new, similar probes can be added with minimal changes elsewhere.

---
