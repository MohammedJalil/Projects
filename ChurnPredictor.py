# -*- coding: utf-8 -*-
"""
AI-Powered Customer Churn Prediction Project
Steps: Load Data → EDA → Preprocess → Train Models → Evaluate → Deploy Demo
"""

# Step 0: Install dependencies (run once)
# !pip install pandas numpy matplotlib seaborn scikit-learn xgboost h2o streamlit

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.metrics import classification_report, roc_auc_score, confusion_matrix
from xgboost import XGBClassifier
import h2o
from h2o.automl import H2OAutoML
import streamlit as st

# ========================
# STEP 1: Load & Explore Data
# ========================
def load_data():
    # Download from Kaggle: https://www.kaggle.com/datasets/blastchar/telco-customer-churn
    df = pd.read_csv('WA_Fn-UseC_-Telco-Customer-Churn.csv')
    
    # Fix TotalCharges (11 empty strings → NaN → fill with median)
    df['TotalCharges'] = pd.to_numeric(df['TotalCharges'], errors='coerce')
    df['TotalCharges'].fillna(df['TotalCharges'].median(), inplace=True)
    
    return df

df = load_data()

# Quick EDA
print("\n=== Data Overview ===")
print(df.head())
print(f"\nMissing Values:\n{df.isnull().sum()}")

# Plot churn distribution
plt.figure(figsize=(8, 4))
sns.countplot(data=df, x='Churn')
plt.title("Churn Distribution (Imbalanced Dataset)")
plt.show()

# ========================
# STEP 2: Preprocess Data
# ========================
def preprocess_data(df):
    # Separate features and target
    X = df.drop(['customerID', 'Churn'], axis=1)
    y = df['Churn'].map({'Yes': 1, 'No': 0})
    
    # Define preprocessing steps
    numeric_features = ['tenure', 'MonthlyCharges', 'TotalCharges']
    categorical_features = ['Contract', 'PaymentMethod', 'InternetService']
    
    preprocessor = ColumnTransformer(
        transformers=[
            ('num', StandardScaler(), numeric_features),
            ('cat', OneHotEncoder(drop='first'), categorical_features)
        ])
    
    X_processed = preprocessor.fit_transform(X)
    return X_processed, y, preprocessor

X, y, preprocessor = preprocess_data(df)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# ========================
# STEP 3: Train Models
# ========================
# Model 1: XGBoost
xgb_model = XGBClassifier(scale_pos_weight=np.sum(y == 0) / np.sum(y == 1))  # Handle class imbalance
xgb_model.fit(X_train, y_train)

# Model 2: AutoML (H2O.ai)
h2o.init()
h2o_df = h2o.H2OFrame(pd.concat([pd.DataFrame(X), y], axis=1))
train, test = h2o_df.split_frame(ratios=[0.8])

aml = H2OAutoML(max_models=5, seed=42)
aml.train(y='Churn', training_frame=train)

# ========================
# STEP 4: Evaluate Models
# ========================
def evaluate_model(model, X_test, y_test, model_name):
    y_pred = model.predict(X_test)
    y_proba = model.predict_proba(X_test)[:, 1] if hasattr(model, 'predict_proba') else None
    
    print(f"\n=== {model_name} Performance ===")
    print(classification_report(y_test, y_pred))
    if y_proba is not None:
        print(f"ROC-AUC: {roc_auc_score(y_test, y_proba):.2f}")
    
    # Plot feature importance (XGBoost only)
    if model_name == "XGBoost":
        plt.figure(figsize=(10, 6))
        feat_imp = pd.Series(model.feature_importances_, 
                            index=preprocessor.get_feature_names_out())
        feat_imp.nlargest(10).plot(kind='barh')
        plt.title("Top 10 Features Driving Churn (XGBoost)")
        plt.show()

evaluate_model(xgb_model, X_test, y_test, "XGBoost")
best_aml_model = aml.leader
evaluate_model(best_aml_model, test, "H2O AutoML")

# ========================
# STEP 5: Streamlit Demo (Optional)
# ========================
def run_streamlit_demo():
    st.title("🔮 Customer Churn Predictor")
    st.write("Predict if a customer will churn using AI")
    
    # User inputs
    tenure = st.slider("Tenure (months)", 0, 72, 12)
    monthly_charges = st.slider("Monthly Charges ($)", 0, 200, 70)
    contract = st.selectbox("Contract", ["Month-to-month", "One year", "Two year"])
    
    if st.button("Predict Churn Risk"):
        # Preprocess input
        input_data = pd.DataFrame({
            'tenure': [tenure],
            'MonthlyCharges': [monthly_charges],
            'TotalCharges': [tenure * monthly_charges],
            'Contract': [contract],
            'PaymentMethod': ["Electronic check"],  # Default value
            'InternetService': ["Fiber optic"]      # Default value
        })
        input_processed = preprocessor.transform(input_data)
        
        # Predict
        churn_prob = xgb_model.predict_proba(input_processed)[0][1]
        st.success(f"Churn Risk: {churn_prob:.0%}")
        st.write("High risk?" if churn_prob > 0.5 else "Low risk")

# Uncomment to run the demo:
# run_streamlit_demo()