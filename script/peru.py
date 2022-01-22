import pandas as pd
import json
import ssl

from urllib.request import Request, urlopen
from datetime import datetime

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

print(last_modified)
print(test_url)

df = pd.read_csv(test_url, delimiter="|")

print(df)