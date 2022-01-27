import pandas as pd

source_url = "https://raw.githubusercontent.com/jmcastagnetto/covid-19-peru-data/main/datos/covid-19-peru-data.csv"

df = pd.read_csv(source_url, delimiter=",")
df = df.rename(columns={"date": "Date", "confirmed": "tot_positive", "total_tests": "tot_tests"})
df = df[df["region"].isna() == True]

df = df[["Date", "tot_positive", "tot_tests"]].reset_index()
df['Date'] = pd.to_datetime(df['Date'], format="%Y-%m-%d")

# df = df.groupby(['Date'], as_index=False).sum()

df["positive"] = 0
df["tests"] = 0


def calc(x):
    df.loc[[x.name], "tests"] = x.tot_tests if x.name == 0 else x.tot_tests - df.iloc[x.name - 1]["tot_tests"]
    df.loc[[x.name], "positive"] = x.tot_positive if x.name == 0 else x.tot_positive - df.iloc[x.name - 1][
        "tot_positive"]


df.apply(calc, axis=1)

df = df.rename(columns={"tests": "Daily change in cumulative total"}).fillna(0)
df = df[df["Daily change in cumulative total"] != 0]

df["Positive rate"] = (
        df["positive"].rolling(7).sum() / df["Daily change in cumulative total"].rolling(7).sum()
).round(3)

df["Positive rate"] = df["Positive rate"].fillna(0)
df["Daily change in cumulative total"] = df["Daily change in cumulative total"].astype(int)

df = df[["Date", "Daily change in cumulative total", "Positive rate"]]

df["Country"] = "Peru"
df["Units"] = "tests performed"
df["Source URL"] = "https://github.com/jmcastagnetto/covid-19-peru-data"
df["Source label"] = "National Institute of Health"
df["Notes"] = pd.NA

df.to_csv("Peru_alt.csv", index=False)
