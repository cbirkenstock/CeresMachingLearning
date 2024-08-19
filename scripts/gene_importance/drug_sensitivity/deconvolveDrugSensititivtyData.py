import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import Lasso
from sklearn.metrics import mean_absolute_error
from sklearn.linear_model import LassoCV
import numpy as np

drug_sensitivity_df = pd.read_csv("given_data/Drug_Sensitivity.csv").dropna()
drug_sensitivity_df["Targets"] = drug_sensitivity_df["Targets"].str.split(",")
drug_sensitivity_df["Targets"] = drug_sensitivity_df["Targets"].apply(
    lambda target_list: [target.replace(" ", "") for target in target_list]
)


target_dummies = (
    pd.get_dummies(drug_sensitivity_df["Targets"].apply(pd.Series).stack())
    .groupby(level=0)
    .sum()
)

drug_sensitivity_df = pd.concat([drug_sensitivity_df, target_dummies], axis=1)
drug_sensitivity_df = drug_sensitivity_df.drop(
    columns=["ID", "Drug Name", "Targets", "Count"]
)


indep_var_array = drug_sensitivity_df.drop(columns=["Z Score"], axis=1).to_numpy()
dep_var_array = drug_sensitivity_df["Z Score"].to_numpy()

ans = np.linalg.lstsq(indep_var_array, dep_var_array, rcond=None)

coefficient_df = pd.DataFrame(
    {
        "geneSymbol": drug_sensitivity_df.drop(columns=["Z Score"], axis=1).columns,
        "Coefficient": ans[0],
    }
)

coefficient_df = coefficient_df[coefficient_df["geneSymbol"] != ""]

k562_df = pd.read_csv("given_data/K562-DHS-v6.noChIP.csv")


merged_df = pd.merge(k562_df, coefficient_df, on="geneSymbol", how="inner")

correlation = merged_df["wgCERES_score"].corr(merged_df["Coefficient"])
print(correlation)

# print(len(merged_df))
# print(merged_df[["wgCERES_score", "Coefficient"]])

# x = drug_sensitivity_df.drop(columns=["Z Score"], axis=1)
# y = drug_sensitivity_df["Z Score"]
# x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=0)

# lr = LinearRegression().fit(x_train, y_train)c

# # Predict and evaluate
# y_pred = lr.predict(x_train)
# mse = mean_squared_error(y_train, y_pred)
# print(f"Mean Squared Error on Test Data: {mse:.2f}")


# # Using LassoCV to find the optimal alpha value with cross-validation
# lasso_cv = LassoCV(alphas=None, cv=10, max_iter=10000)
# lasso_cv.fit(x_train, y_train)

# alpha_value = lasso_cv.alpha_
# lasso = Lasso(alpha=alpha_value, max_iter=10000)
# lasso.fit(x_train, y_train)

# y_pred_lasso = lasso.predict(x_test)
# mse_lasso = mean_squared_error(y_test, y_pred_lasso)
# print(f"Mean Squared Error with Lasso Regression: {mse_lasso}")
