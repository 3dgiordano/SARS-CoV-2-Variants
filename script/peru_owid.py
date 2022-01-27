import pandas as pd
from cowidev.testing import CountryTestBase
from cowidev.utils.web import request_json


class Peru(CountryTestBase):
    location = "Peru"
    units = "tests performed"
    source_label = "National Institute of Health"
    source_url = "https://datos.ins.gob.pe/api/3/action/package_show?id=dataset-de-pruebas-moleculares-del-instituto-nacional-de-salud-ins"
    source_url_ref = "https://datos.ins.gob.pe/dataset/dataset-de-pruebas-moleculares-del-instituto-nacional-de-salud-ins"

    def read(self) -> pd.DataFrame:
        json_dict = request_json(self.source_url)
        resources = json_dict["result"]["resources"]
        last_modified = max(datetime.fromisoformat(node["last_modified"]) for node in resources)
        test_url = [obj["url"] for obj in resources if datetime.fromisoformat(obj["last_modified"]) == last_modified][0]
        return pd.read_csv(test_url, delimiter="|", low_memory=False)

    @staticmethod
    def pipe_normalize(df: pd.DataFrame) -> pd.DataFrame:
        df = df.dropna()
        df = df.rename(columns={"FECHA_MUESTRA": "Date", "RESULTADO": "Result"})
        df = df[["Date", "Result"]]
        df = df[df["Date"] >= 20200101]
        df['Date'] = pd.to_datetime(df['Date'], format="%Y%m%d")

        df["positive"] = df["Result"].str.contains("POS").astype(int)
        df["negative"] = df["Result"].str.contains("NEG").astype(int)

        df = df[["Date", "positive", "negative"]]

        df = df.groupby(['Date'], as_index=False).sum()
        return df

    @staticmethod
    def pipe_metrics(df: pd.DataFrame) -> pd.DataFrame:
        df["Daily change in cumulative total"] = df["positive"] + df["negative"]
        df["Daily change in cumulative total"] = df["Daily change in cumulative total"].astype(int)

        df = df[df["Daily change in cumulative total"] != 0]

        df["Positive rate"] = (
                df["positive"].rolling(7).sum() / df["Daily change in cumulative total"].rolling(7).sum()
        ).round(3)

        df["Positive rate"] = df["Positive rate"].fillna(0)

        df = df[["Date", "Daily change in cumulative total", "Positive rate"]]
        return df

    def pipeline(self, df: pd.DataFrame) -> pd.DataFrame:
        return df.pipe(self.pipe_normalize).pipe(self.pipe_metrics).pipe(self.pipe_metadata)

    def export(self):
        df = self.read().pipe(self.pipeline)
        self.export_datafile(df)


def main():
    Peru().export()


if __name__ == "__main__":
    main()
