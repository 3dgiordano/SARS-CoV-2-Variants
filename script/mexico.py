import datetime

import pandas as pd

yesterday = (datetime.date.today() - datetime.timedelta(days=1)).strftime("%Y%m%d")

base_url = "https://datos.covid-19.conacyt.mx/Downloads/Files/Casos_Diarios_Estado_Nacional_"
confirmed_url = f"{base_url}Confirmados_{yesterday}.csv"
negatives_url = f"{base_url}Negativos_{yesterday}.csv"


def url_melt(url, name):
    df_melt = pd.read_csv(url).melt(id_vars=["cve_ent", "poblacion", "nombre"], var_name="Date", value_name=name)
    df_melt = df_melt[df_melt["nombre"] == "Nacional"]
    df_melt = df_melt.drop(["cve_ent", "poblacion", "nombre"], axis=1)
    df_melt['Date'] = pd.to_datetime(df_melt['Date'], format="%d-%m-%Y")
    return df_melt


df_confirmed = url_melt(confirmed_url, "positive")

df_negatives = url_melt(negatives_url, "negative")

df = pd.merge(df_confirmed, df_negatives, on=["Date"], how="right").fillna(0).sort_values("Date")

df["Daily change in cumulative total"] = df["positive"] + df["negative"]
df["Daily change in cumulative total"] = pd.to_numeric(df["Daily change in cumulative total"], downcast='integer')

df = df[df["Daily change in cumulative total"] != 0]

df["Positive rate"] = (
        df["positive"].rolling(7).sum() / df["Daily change in cumulative total"].rolling(7).sum()
).round(3)

df["Positive rate"] = df["Positive rate"].fillna(0)

df = df[["Date", "Daily change in cumulative total", "Positive rate"]]

df["Country"] = "Mexico"
df["Units"] = "people tested"
df["Source URL"] = "https://datos.covid-19.conacyt.mx/#DownZCSV"
df["Source label"] = "Health Secretary"
df["Notes"] = pd.NA

df.to_csv("Mexico_wip.csv", index=False)
