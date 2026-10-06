# Exact-arithmetic verification

Checks for [*Pinching cones for positive isotropic curvature in dimensions seven and eight*](https://arxiv.org/abs/2608.26598).

[View the notebook with saved outputs](https://github.com/JaehoCho43/pic8-verification/blob/main/verification.ipynb). No installation is needed to view it on GitHub.

## Run the Python script

Requires Python 3 and SymPy. The recorded run used Python 3.14.7 and SymPy 1.14.0.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install sympy==1.14.0
python verify_exact.py
```

On Windows, activate the environment with `.venv\Scripts\activate` instead.

Successful execution ends with:

```text
ALL 178 EXACT CHECKS PASSED. No geometric proof claim is made.
```

## Run the notebook

In the same environment:

```bash
python -m pip install jupyterlab ipykernel
python -m ipykernel install --user --name pic8-verification --display-name "PIC verification"
jupyter lab verification.ipynb
```

Select the **PIC verification** kernel, then **Restart Kernel / Run All**. Saved outputs are already included in the notebook. Each cell has a short description of the relevant part of the paper and the checking method.

`verify_exact.py` and the mathematical code cells in `verification.ipynb` perform the same checks. The notebook also has a cell that records the Python/SymPy versions.

These are checks of selected algebraic calculations, not a verification of the entire geometric proof. No manuscript files or external datasets are required or included.
