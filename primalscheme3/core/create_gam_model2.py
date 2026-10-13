# train_surrogate.py
import joblib
import pandas as pd
import xgboost as xgb
from sklearn.compose import ColumnTransformer
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


def load_data(filepath):
    df = pd.read_csv(filepath)
    X = df.drop(columns=["predicted_delta_ct"])
    y = df["predicted_delta_ct"]
    return train_test_split(X, y, test_size=0.2, random_state=42)


def build_pipeline():
    categorical_cols = ["mut_base", "dinu2_gc_rich"]
    numeric_cols = ["p3_dg", "prdg", "idx_prop", "mm_dg2"]

    preprocessor = ColumnTransformer(
        transformers=[
            ("num", "passthrough", numeric_cols),
            ("cat", OneHotEncoder(handle_unknown="ignore"), categorical_cols),
        ]
    )

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
                    n_jobs=-1,
                ),
            ),
        ]
    )
    return pipeline


def train_and_evaluate(csv_filepath, model_save_path):
    X_train, X_test, y_train, y_test = load_data(csv_filepath)

    pipeline = build_pipeline()
    print("Training XGBoost surrogate...")
    pipeline.fit(X_train, y_train)

    predictions = pipeline.predict(X_test)
    r2 = r2_score(y_test, predictions)
    rmse = mean_squared_error(y_test, predictions) ** 0.5

    print("\nSurrogate Performance:")
    print(f"  R² = {r2:.4f} (target > 0.98)")
    print(f"  RMSE = {rmse:.4f} ΔCt")

    joblib.dump(pipeline, model_save_path)
    print(f"Saved to {model_save_path}")


if __name__ == "__main__":
    train_and_evaluate(
        "/Users/kentcg/qpcr-single-mm/dinucleotide/q5_gam_surrogate_training_data.csv",
        "q5_gam_surrogate_training_data.joblib",
    )
