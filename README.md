# Phylogenetic Tree Simulation and Analysis

**Author:** Mattia Colbertaldo  
**Date:** February 2024

---

[![Work in Progress](https://img.shields.io/badge/Work%20in%20Progress-Yes-yellow.svg)](https://shields.io/)


## Description

This repository contains code and data for a thesis project focused on simulating phylogenetic trees under different models and inferring their parameters using machine learning (ML) techniques.

The main purpose of this thesis work is to substitute the slow Maximum Likelihood Estimation (MLE) methods traditionally used to infer parameters with ML methods. Additionally, the goal is to create a model predictor that, given a phylogenetic tree, can automatically guess the model and infer its parameters.

### Workflow

1. **Simulations:** 
    - The `simulations.r` file is used to simulate phylogenetic trees under different models using functions from the diversitree package. The trees are then saved to file for further analysis.

2. **Parameter Inference:**
    - In the `ranges.r` file, ranges of parameters for simulating trees are obtained using inference on phylogenetic trees from real-world data.
    - `CDV_full_tree.py` encodes the trees with the CDV representation, preparing them for ML analysis.
    - In the `AllModels.ipynb` notebook, encoded trees are read, data is managed, and they are input into a Convolutional Neural Network (CNN) created to train it to infer the parameters of the model. The code is universal for all models in the diversitree package.

3. **Model Predictor:**
    - The `ModelPredictor.ipynb` notebook trains different Neural Networks to try to infer the model given a tree.

---

## Models

We explore various models for simulating phylogenetic trees, including:
- BD
- BiSSE
- MuSSE
- QuaSSE
- GeoSSE
- BiSSEness
- ClaSSE

Additionally, we compare simulations with the Constant Birth-Death model.

---

## Usage

Feel free to explore the code and data provided in this repository. For detailed instructions on running simulations, inferring parameters, and training model predictors, refer to the respective files and notebooks.

---

## References

For further reading and understanding of the models, methods, and ML techniques used, please refer to the cited references.

---

If you have any questions or suggestions, don't hesitate to reach out. Happy exploring!

