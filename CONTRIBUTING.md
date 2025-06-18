# Contributing to Sandy

First off, thanks for taking the time to contribute to **Sandy**! 🎉  
This document outlines the guidelines for contributing to this project.

## 🙌 Who Can Contribute?

**Everyone is welcome to contribute.** Whether it's fixing bugs, improving documentation, or proposing new features, your help is appreciated.

---

## 🚀 Branching & Versioning Strategy

- The **main branch is `develop`**.
- All new version developments are done in branches named **`vX.Y`**, where `X.Y` indicates the upcoming version (e.g., `v1.1`).
- All **features or fixes are branched from `vX.Y`** and submitted via **Pull Requests** (PRs).
- **Tags and full releases to PyPI** are run via action workflows.
- **GitHub releases are added manually**.

**Current active development branch:** `v1.1`

---

## 🧱 Code Structure and Design

- **Keep changes minimal** and **object-oriented**.
- Follow existing architectural patterns—prefer adding to existing modules unless there's a compelling reason to create new ones.

---

## 📚 Docstrings & Testing

Each **class, method, and function must**:
- Use **NumPy-style docstrings**.
- Include **working examples** inside the docstring.
- Examples will be automatically run as part of the test suite via `pytest`.

### Example Docstring

```python
def square(x):
    """
    Return the square of a number.

    Parameters
    ----------
    x : int or float
        The number to square.

    Returns
    -------
    int or float
        The square of `x`.

    Examples
    --------
    >>> square(2)
    4
    >>> square(-3)
    9
    """
    return x * x
```

# 🧪 Running Tests

To run all tests, including docstring examples:

```bash
pytest --doctest-modules
```

Make sure all tests pass before submitting a pull request.

## 📝 Submitting a Pull Request

1. **Fork the repository** and create a branch from `vX.Y`.
2. **Make your changes**, ensuring each addition includes:
   - NumPy-style docstrings
   - Working examples in the docstrings
   - Associated tests
3. **Run the test suite** to make sure all tests pass:

   ```bash
   pytest --doctest-modules
   ``

4. **Submit your pull request** with a short but descriptive title.

5. **If you bump the version**, use the appropriate format:
   - X.Y-betaN for changes on vX.Y branches (beta releases)
   - X.Y for changes merged into develop (official release)

---

Thank you for helping improve Sandy! 🏖️