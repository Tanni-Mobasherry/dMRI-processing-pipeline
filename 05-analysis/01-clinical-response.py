import pandas as pd

# Load clinical data
df = pd.read_excel("clinical-output.xlsx")

# Keep participants with both BL and FU DASS scores
df = df.dropna(subset=["ID", "BL_MRIdepDASS", "FU_MRIdepDASS"]).copy()

# Calculate DASS reduction
df["DASS_Reduction"] = (
    (df["BL_MRIdepDASS"] - df["FU_MRIdepDASS"])
    / df["BL_MRIdepDASS"]
)

# Responder = at least 50% reduction
df["Responder_Status"] = df["DASS_Reduction"].apply(
    lambda x: "Responder" if x >= 0.50 else "Non-responder"
)

# Keep only what we need
result = df[
    [
        "ID",
        "BL_MRIdepDASS",
        "FU_MRIdepDASS",
        "DASS_Reduction",
        "Responder_Status",
    ]
]

# Save
result.to_csv("clinical_response.csv", index=False)

# Show result
print(result.to_string(index=False))

print("\nResponders:",
      (result["Responder_Status"] == "Responder").sum())

print("Non-responders:",
      (result["Responder_Status"] == "Non-responder").sum())
