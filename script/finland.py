import pandas as pd

base_url = "https://sampo.thl.fi/pivot/prod/en/epirapo/covid19case/fact_epirapo_covid19case.csv"
test_url = f"{base_url}?row=dateweek20200101-509093L&column=measure-444833.445356.492118.&&fo=1"

df = pd.read_csv(test_url, delimiter=";")
df['Time'] = pd.to_datetime(df['Time'], format="%Y-%m-%d")

df = df.pivot(index='Time', columns='Measure', values='val').fillna(0).reset_index()
df = df.rename(
    columns={"Time": "Date", "Number of cases": "positive", "Number of tests": "Daily change in cumulative total"})
df = df.drop("Number of deaths", axis=1)

df["Daily change in cumulative total"] = pd.to_numeric(df["Daily change in cumulative total"], downcast='integer')
df = df[df["Daily change in cumulative total"] != 0]

df["Positive rate"] = (
        df["positive"].rolling(7).sum() / df["Daily change in cumulative total"].rolling(7).sum()
).round(3)

df["Positive rate"] = df["Positive rate"].fillna(0)

df = df[["Date", "Daily change in cumulative total", "Positive rate"]]

df["Country"] = "Finland"
df["Units"] = "tests performed"
df["Source URL"] = "https://sampo.thl.fi/pivot/prod/en/epirapo/covid19case/fact_epirapo_covid19case"
df["Source label"] = "Finnish Institute for Health and Welfare"
df["Notes"] = pd.NA

print(df)

df.to_csv("Finland_wip.csv", index=False)
