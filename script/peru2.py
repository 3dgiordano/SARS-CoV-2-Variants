import json
import ssl
from datetime import datetime
from urllib.request import Request, urlopen

import pandas as pd

source_url = "https://www.datosabiertos.gob.pe/api/3/action/package_show?id=dataset-de-pruebas-moleculares-del-instituto-nacional-de-salud-para-covid-19-ins"

context = ssl._create_unverified_context()

ssl._create_default_https_context = ssl._create_unverified_context

body_json = urlopen(
    Request(source_url, headers={'User-Agent': 'Mozilla/5.0'}), context=context
).read().decode('UTF-8')
json_dict = json.loads(body_json)
resources = json_dict["result"][0]["resources"]

last_modified = max(datetime.strptime(node["revision_timestamp"].split(",")[1].strip(), "%m/%d/%Y - %H:%M") for node in resources)
test_url = [obj["url"] for obj in resources if
            datetime.strptime(obj["revision_timestamp"].split(",")[1].strip(), "%m/%d/%Y - %H:%M") == last_modified][0]

print(test_url)

df = pd.read_csv(test_url, delimiter="|")

print(df)

