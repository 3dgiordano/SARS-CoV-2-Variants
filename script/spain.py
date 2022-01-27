import datetime

import pandas as pd

yesterday = (datetime.date.today() - datetime.timedelta(days=1)).strftime("%d%m%Y")

base_url = "https://www.sanidad.gob.es/profesionales/saludPublica/ccayes/alertasActual/nCov"
test_url = f"{base_url}/documentos/Datos_Pruebas_Realizadas_Historico_{yesterday}.csv"

df = pd.read_csv(test_url, encoding='cp1252', delimiter=";")
df = df.drop(["PROVINCIA"], axis=1)
df = df.rename(columns={"FECHA_PRUEBA": "Date"})
df['Date'] = pd.to_datetime(df['Date'], format="%d%b%Y")

df = df.groupby(['Date'])[['N_ANT_POSITIVOS', 'N_PCR_POSITIVOS', 'N_ANT', 'N_PCR']].apply(sum).reset_index()

df["positive"] = df["N_ANT_POSITIVOS"] + df["N_PCR_POSITIVOS"]
df["Daily change in cumulative total"] = df["N_ANT"] + df["N_PCR"]
df["Daily change in cumulative total"] = pd.to_numeric(df["Daily change in cumulative total"], downcast='integer')

df = df[df["Daily change in cumulative total"] != 0]

df["Positive rate"] = (
        df["positive"].rolling(7).sum() / df["Daily change in cumulative total"].rolling(7).sum()
).round(3)

df["Positive rate"] = df["Positive rate"].fillna(0)

df = df[["Date", "Daily change in cumulative total", "Positive rate"]]

df["Country"] = "Spain"
df["Units"] = "tests performed"
df["Source URL"] = f"{base_url}/pruebasRealizadas.htm"
df["Source label"] = "Ministry of Health"
df["Notes"] = pd.NA

df.to_csv("Spain_wip.csv", index=False)
