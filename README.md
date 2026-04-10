A Hybrid Transformer and Gradient Boosting Approach for BGC Detection and Functional Genomic Insights

## Overview

This project implements a **hybrid deep learning + machine learning pipeline** for identifying Biosynthetic Gene Clusters (BGCs) and classifying their types using transformer-based feature extraction followed by XGBoost classification, along with **SHAP-based biological interpretation** to identify important k-mers and uncover biologically meaningful genomic patterns.

---

## Workflow Summary

1. 6-mer Feature Extraction
2. Dataset Preparation (BGC vs Non-BGC)
3. Grid Search Hyperparameter Optimization
4. General Transformer Training (Binary Classification)
5. Class-Specific Transformer Training (BGC Types)
6. Feature Extraction & Fusion
7. XGBoost Training + Cross Validation
8. Test Evaluation
9. SHAP-based Interpretation

---

## 1. k-mer Feature Extraction

* DNA sequences converted into k-mer frequency vectors
* Example: k = 6 → 1024 features
* Output: numerical feature matrix

---

## 2. Dataset Preparation

* Input CSV:

  * Feature columns → k-mer features
  * Last column → Label

### Label Structure

* Class 0 → Non-BGC
* Class 1..N → BGC types

### Binary Conversion

* Non-BGC → 0
* BGC → 1

---

## 3. General Transformer (BGC vs Non-BGC)

### Purpose

* Learn global genomic patterns
* Separate BGC vs Non-BGC


---

## 4. Class-Specific Transformers

### Purpose

* One transformer per BGC class
* Learn class-specific patterns

---

## 5. Feature Extraction & Fusion

### Extracted Features

* General transformer features
* Class-specific transformer features

### Final Feature Vector

Final Features = [General Features + All Class Features]

### Advantage

* Captures:

  * Global genomic patterns
  * Class-specific signatures

---

## 6. XGBoost Classification

### Input

* Combined transformer features

### Model Parameters

* n_estimators = 400
* max_depth = 6
* learning_rate = 0.05
* subsample = 0.8
* colsample_bytree = 0.8
* objective = multi:softprob

### Training Strategy

* Stratified 5-Fold Cross Validation

### Metrics

* Accuracy
* Precision (weighted)
* Recall (weighted)
* F1-score
* AUROC

---

## 7. Grid Search Optimization

### Purpose

* Tune Transformer , XGBoost hyperparameters

### Outputs

* Best hyperparameters
* Grid search results
* Optimized model

---

## 8. Test Evaluation

* Final model trained on Train + Validation
* Evaluated on independent Test set

### Outputs

* Test metrics
* Final trained model

---

## 9. SHAP-Based Biological Interpretation

### Purpose

* Interpret model predictions
* Identify important k-mers

### Biological Insights

* Detect enriched sequence patterns
* Link k-mers to functional regions
* Understand genomic signatures of BGCs

### Outputs

* SHAP values
* Feature importance plots
* Top contributing k-mers

---

## Directory Structure

```
project/
│── data/
│   ├── raw_sequences.fasta
│   ├── labels.csv
│   ├── final_output.csv
│
│── models/
│   ├── PRETRAINED_TRANSFORMERS_BGC/
│   │   ├── bgc_general_transformer.keras
│   │   ├── best_bgc_general.keras
│   │   ├── bgc_class_transformer_<class>.keras
│   │   ├── best_bgc_class_<class>.keras
│   │   ├── scaler.pkl
│   │   ├── label_encoder.pkl
│   │
│   ├── xgboost_final_model.pkl
│   
│
│── scripts/
│   ├── kmer_generation.py
│   ├── non_bgc_extraction.py
│   ├── grid_search.py
│   ├── bgc_classification.py
│   ├── shap_analysis.py
│
│── results/
│   ├── TRAINING_OUTPUT/
│   │   ├── X_train_combined.npy
│   │   ├── X_val_combined.npy
│   │   ├── X_test_combined.npy
│   │   ├── y_train.npy
│   │   ├── y_val.npy
│   │   ├── y_test.npy
│   │   ├── train_gen_features.npy
│   │   ├── test_gen_features.npy
│   │   ├── training_metadata.csv
│   │   ├── early_stopping_stats.csv
│   │   ├── class_info.csv
│   │   ├── general_transformer_history.csv
│   │   ├── class_<class>_history.csv
│   │
│   │
│   ├── shap_results/
│   │   ├── shap_values.npy
│   │   ├── shap_summary_plot.png
│   │   ├── shap_feature_importance.csv
│   │   ├── top_kmers.txt

```

---

## Key Highlights

* Hybrid Transformer + XGBoost pipeline
* Memory-optimized architecture
* Early stopping for stability
* Class-specific learning
* Feature fusion strategy
* Grid search optimization
* SHAP-based biological interpretation

---


