import joblib
import pandas as pd
import xgboost as xgb
from sklearn.compose import ColumnTransformer
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


def load_data(filepath):
    """Loads the synthetic data and separates features from the target."""
    print("Loading synthetic data from R...")
    df = pd.read_csv(filepath)

    # X contains everything EXCEPT the target column
    X = df.drop(columns=["predicted_delta_ct"])
    # y is ONLY the target column
    y = df["predicted_delta_ct"]

    return train_test_split(X, y, test_size=0.2, random_state=42)


def build_xgboost_pipeline():
    """Builds a pipeline that automatically encodes text and runs XGBoost."""
    # 1. Define which columns are text and which are numbers
    categorical_cols = ["mut_base", "dinu1", "dinu2"]
    numeric_cols = ["dg", "ddg", "p3_dg", "idx", "idx_prop"]

    # 2. Create the Preprocessor (Translates text like 'TC/AC' into 1s and 0s)
    preprocessor = ColumnTransformer(
        transformers=[
            ("num", "passthrough", numeric_cols),
            ("cat", OneHotEncoder(handle_unknown="ignore"), categorical_cols),
        ]
    )

    # 3. Bundle the Preprocessor and the XGBoost Regressor together
    # These parameters (n_estimators=500, max_depth=6) make the model
    # highly flexible so it can perfectly trace the GAM's non-linear curves.
    pipeline = Pipeline(
        steps=[
            ("preprocessor", preprocessor),
            (
                "regressor",
                xgb.XGBRegressor(
                    n_estimators=500,
                    max_depth=6,
                    learning_rate=0.05,
                    random_state=42,
                    n_jobs=-1,  # Uses all available CPU cores for faster training
                ),
            ),
        ]
    )

    return pipeline


def train_and_evaluate(csv_filepath, model_save_path="primer_gam_surrogate.joblib"):
    """Main function to load data, train the model, and save it to disk."""
    # 1. Get the data
    X_train, X_test, y_train, y_test = load_data(csv_filepath)

    # 2. Build the model pipeline
    model_pipeline = build_xgboost_pipeline()

    # 3. Train it!
    print("Training the XGBoost surrogate model... (this may take a minute)")
    model_pipeline.fit(X_train, y_train)

    # 4. Evaluate its accuracy against the R model's predictions
    predictions = model_pipeline.predict(X_test)
    r2 = r2_score(y_test, predictions)
    rmse = mean_squared_error(y_test, predictions) ** 0.5

    print("\n--- Model Performance ---")
    print(f"R-squared Accuracy: {r2:.4f} (Aiming for > 0.98)")
    print(f"Average Error (RMSE): {rmse:.4f} Delta Ct cycles")
    print("-------------------------\n")

    # 5. Save the trained pipeline
    joblib.dump(model_pipeline, model_save_path)
    print(f"Success! Production model saved to: {model_save_path}")


# ==========================================
# Run the training process!
# ==========================================
if __name__ == "__main__":
    # Ensure this matches the name of the CSV you generated in R
    train_and_evaluate(
        "/Users/kentcg/qpcr-single-mm/dinucleotide/q5_gam_surrogate_training_data.csv",
        "/Users/kentcg/qpcr-single-mm/dinucleotide/q5_gam_surrogate_training_data.joblib",
    )
