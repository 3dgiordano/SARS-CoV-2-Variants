import pandas as pd

from cowidev.testing import CountryTestBase


class Finland(CountryTestBase):
    location = "Finland"
    units = "tests performed"
    base_url = "https://sampo.thl.fi/pivot/prod/en/epirapo/covid19case/fact_epirapo_covid19case.csv"
    source_url = f"{base_url}?row=dateweek20200101-509093L&column=measure-444833.445356.492118.&&fo=1"
    source_url_ref = "https://sampo.thl.fi/pivot/prod/en/epirapo/covid19case/fact_epirapo_covid19case"
    source_label = "Finnish Department of Health and Welfare"

    def read(self) -> pd.DataFrame:
        df = pd.read_csv(source_url, delimiter=";")
        df['Time'] = pd.to_datetime(df['Time'], format="%Y-%m-%d")
        return df

    def pipe_metrics(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df.pivot(index='Time', columns='Measure', values='val').fillna(0).reset_index()
        df = df.rename(
            columns={"Time": "Date", "Number of cases": "positive",
                     "Number of tests": "Daily change in cumulative total"})
        df = df.drop("Number of deaths", axis=1)

        df["Daily change in cumulative total"] = pd.to_numeric(df["Daily change in cumulative total"],
                                                               downcast='integer')
        df = df[df["Daily change in cumulative total"] != 0]

        df["Positive rate"] = (
                df["positive"].rolling(7).sum() / df["Daily change in cumulative total"].rolling(7).sum()
        ).round(3)

        df["Positive rate"] = df["Positive rate"].fillna(0)
        return df

    def pipe_filter_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        return df[["Date", "Daily change in cumulative total", "Positive rate"]]

    def pipeline(self, df: pd.DataFrame) -> pd.DataFrame:
        return df.pipe(self.pipe_metrics).pipe(self.pipe_filter_columns).pipe(self.pipe_metadata)

    def export(self):
        df = self.read().pipe(self.pipeline)
        self.export_datafile(df)


def main():
    Finland().export()


if __name__ == "__main__":
    main()
