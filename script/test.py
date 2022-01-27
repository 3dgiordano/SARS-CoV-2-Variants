import datetime

import pandas as pd
from cowidev.testing import CountryTestBase


class Mexico(CountryTestBase):
    location = "Mexico"
    units = "people tested"
    source_label = "Health Secretary"
    source_url = "https://datos.covid-19.conacyt.mx/#DownZCSV"
    source_url_ref = source_url

    def url_melt(self, url: str, name: str) -> pd.DataFrame:
        df_melt = pd.read_csv(url).melt(id_vars=["cve_ent", "poblacion", "nombre"], var_name="Date", value_name=name)
        df_melt = df_melt[df_melt["nombre"] == "Nacional"]
        df_melt = df_melt.drop(["cve_ent", "poblacion", "nombre"], axis=1)
        df_melt['Date'] = pd.to_datetime(df_melt['Date'], format="%d-%m-%Y")
        return df_melt

    def read(self) -> pd.DataFrame:
        yesterday = (datetime.date.today() - datetime.timedelta(days=1)).strftime("%Y%m%d")
        base_url = "https://datos.covid-19.conacyt.mx/Downloads/Files/Casos_Diarios_Estado_Nacional_"
        return pd.merge(
            self.url_melt(f"{base_url}Confirmados_{yesterday}.csv", "positive"),
            self.url_melt(f"{base_url}Negativos_{yesterday}.csv", "negative"),
            on=["Date"], how="right").fillna(0).sort_values("Date")

    def pipe_daily_change_in_cum_total(self, df: pd.DataFrame) -> pd.DataFrame:
        df["Daily change in cumulative total"] = df["positive"] + df["negative"]
        df["Daily change in cumulative total"] = pd.to_numeric(df["Daily change in cumulative total"],
                                                               downcast='integer')
        df = df[df["Daily change in cumulative total"] != 0]
        return df

    def pipe_positive_rate(self, df: pd.DataFrame) -> pd.DataFrame:
        df["Positive rate"] = (
                df["positive"].rolling(7).sum() / df["Daily change in cumulative total"].rolling(7).sum()
        ).round(3)
        df["Positive rate"] = df["Positive rate"].fillna(0)
        return df

    def pipe_filter_rows(self, df: pd.DataFrame) -> pd.DataFrame:
        return df[["Date", "Daily change in cumulative total", "Positive rate"]]

    def pipeline(self, df: pd.DataFrame) -> pd.DataFrame:
        return (
            df.pipe(self.pipe_daily_change_in_cum_total).pipe(self.pipe_positive_rate).pipe(self.pipe_filter_rows).pipe(
                self.pipe_metadata)
        )

    def export(self):
        df = self.read().pipe(self.pipeline)
        self.export_datafile(df)


def main():
    Mexico().export()


if __name__ == "__main__":
    main()
