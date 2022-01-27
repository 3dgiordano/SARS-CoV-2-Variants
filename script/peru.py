import json
import ssl
from datetime import datetime
from urllib.request import Request, urlopen

import pandas as pd

source_url = "https://datos.ins.gob.pe/api/3/action/package_show?id=dataset-de-pruebas-moleculares-del-instituto-nacional-de-salud-ins"

context = ssl._create_unverified_context()

ssl._create_default_https_context = ssl._create_unverified_context

body_json = urlopen(
   Request(source_url, headers={'User-Agent': 'Mozilla/5.0'}), context=context
).read().decode('UTF-8')
json_dict = json.loads(body_json)
resources = json_dict["result"]["resources"]

last_modified = max(datetime.fromisoformat(node["last_modified"]) for node in resources)
test_url = [obj["url"] for obj in resources if datetime.fromisoformat(obj["last_modified"]) == last_modified][0]

# test_url = "C:\\Users\\David\\Downloads\\pm14enero2022.zip"

df = pd.read_csv(test_url, delimiter="|", low_memory=False)
df = df.dropna()
df = df.rename(columns={"FECHA_MUESTRA": "Date", "RESULTADO": "Result"})
df = df[["Date", "Result"]]
df = df[df["Date"] >= 20200101]
df['Date'] = pd.to_datetime(df['Date'], format="%Y%m%d")

df["positive"] = df["Result"].str.contains("POS").astype(int)
df["negative"] = df["Result"].str.contains("NEG").astype(int)

df = df[["Date", "positive", "negative"]]

df = df.groupby(['Date'], as_index=False).sum()

df["Daily change in cumulative total"] = df["positive"] + df["negative"]
df["Daily change in cumulative total"] = df["Daily change in cumulative total"].astype(int)

df = df[df["Daily change in cumulative total"] != 0]

df["Positive rate"] = (
        df["positive"].rolling(7).sum() / df["Daily change in cumulative total"].rolling(7).sum()
).round(3)

df["Positive rate"] = df["Positive rate"].fillna(0)

df = df[["Date", "Daily change in cumulative total", "Positive rate"]]

df["Country"] = "Peru"
df["Units"] = "tests performed"
df["Source URL"] = "https://datos.ins.gob.pe/dataset/dataset-de-pruebas-moleculares-del-instituto-nacional-de-salud-ins"
df["Source label"] = "National Institute of Health"
df["Notes"] = pd.NA

df.to_csv("Peru_wip.csv", index=False)
