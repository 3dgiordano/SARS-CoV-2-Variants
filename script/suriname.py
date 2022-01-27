import time
from datetime import datetime, timedelta
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import pandas as pd

# https://covid-19.sr/
# Use Wayback machine to get the historical information
source_url = "http://web.archive.org/web/timemap/link/https://covid-19.sr/"

body = urlopen(
    Request(source_url, headers={'User-Agent': 'Mozilla/5.0'})
).read().decode('UTF-8')

sources_url = [
    x.replace("<", "").replace(">", "").replace(' rel="', '').replace('",', '').replace(' type="', '').replace(
        ' from="', '').replace(' datetime="', '').split(";") for x in body.splitlines()]

sources_url = [{'url': x[0], 'date': x[2]} for x in sources_url if "memento" in x[1]]

sources_url.append(
    {'url': "https://covid-19.sr/",
     'date': (datetime.now() - timedelta(days=1)).strftime('%a, %d %b %Y 23:59:59 GMT')}
)

# sources_url = sources_url[-10:]

df_urls = pd.DataFrame.from_records(sources_url)
df_urls['date'] = pd.to_datetime(df_urls['date'])
df_urls['d'] = df_urls['date'].dt.date

max_dates = [str(d) for d in df_urls.groupby(['d'])["date"].max().tolist()]

df_urls["tests"] = 0
df_urls["positive"] = 0

from_date = pd.Timestamp(2020, 8, 3)

dates_excludes = ["2020-10-05 19:46:00+00:00", "2021-01-11 20:22:02+00:00", "2021-03-22 19:01:13+00:00",
                  "2021-04-12 22:41:50+00:00", "2021-04-19 14:05:56+00:00", "2021-11-18 17:48:44+00:00",
                  "2021-05-29 18:12:57+00:00", "2021-09-13 18:52:43+00:00", "2021-08-30 17:46:05+00:00",
                  "2020-10-14 17:37:45+00:00", "2020-10-21 22:18:54+00:00", "2021-03-21 22:26:49+00:00",
                  "2021-10-04 21:09:32+00:00", "2021-10-15 22:18:03+00:00", "2021-10-23 19:48:51+00:00",
                  "2021-11-15 19:03:24+00:00", "2021-11-27 14:28:54+00:00", "2021-11-29 18:53:27+00:00",
                  "2021-11-30 17:10:31+00:00", "2021-12-01 16:34:11+00:00", "2022-01-02 00:48:20+00:00",
                  "2022-01-16 13:24:42+00:00", "2022-01-19 03:55:46+00:00", "2021-11-01 20:54:10+00:00",
                  "2021-10-29 18:24:52+00:00", "2021-10-24 04:54:03+00:00", "2021-09-22 20:36:42+00:00",
                  "2020-12-18 17:59:27+00:00"
                  ]


def iterate(x):
    if str(x.date) in max_dates and str(x.date) not in dates_excludes and x.date.value >= from_date.value:
        tries = 0
        body_url = ""
        redirected = False
        while True:
            try:
                req = urlopen(
                    Request(x.url, headers={'User-Agent': 'Mozilla/5.0'})
                )
                # Don't allow redirect
                if req.url == x.url:
                    body_url = req.read().decode('UTF-8')
                else:
                    redirected = True
                break
            except HTTPError as err:
                print("Retry " + x.url)
                tries += 1
                if tries > 3:
                    break

        tests = 0
        if "Totaal Testen" in body_url:
            tests = int(body_url.split('Totaal Testen')[0].split('data-counter-value="')[-1].split('"')[0])
        test_neg = 0
        if "Totaal negatieve" in body_url:
            test_neg = int(body_url.split('Totaal negatieve')[0].split('data-counter-value="')[-1].split('"')[0])

        positive = tests - test_neg
        if tests > 0 and positive <= tests:
            print(x.url)
            print(str(x.date) + " Test:" + str(tests) + " Positive:" + str(positive) + " Neg:" + str(test_neg))
            df_urls.loc[[x.name], "tests"] = tests
            df_urls.loc[[x.name], "positive"] = positive
        else:
            if not redirected:
                print(
                    "Analyze:" + str(x.url) + " " + str(x.date) + " Test:" + str(tests) + " Positive:" + str(positive))
            else:
                print(
                    "Redirect:" + str(x.url) + " " + str(x.date))


df_urls.apply(iterate, axis=1)

# Clean
df = df_urls[df_urls["tests"] != 0].reset_index()

df["delete"] = 0


def clean_duplicates(x):
    if x.name > 0:
        if df.iloc[x.name - 1]["tests"] == x.tests and df.iloc[x.name - 1]["positive"] == x.positive:
            df.loc[[x.name], "delete"] = 1


df.apply(clean_duplicates, axis=1)
df = df[df["delete"] == 0]

df = df[["date", "tests", "positive"]]

df = df.rename({"date": "Date", "tests": "Daily change in cumulative total"}, axis=1)

df["Positive rate"] = (
        df["positive"].rolling(7).sum() / df["Daily change in cumulative total"].rolling(7).sum()
).round(3)

df["Positive rate"] = df["Positive rate"].fillna(0)

df = df[["Date", "Daily change in cumulative total", "Positive rate"]]

df["Country"] = "Suriname"
df["Units"] = "tests performed"
df["Source URL"] = "https://covid-19.sr/"
df["Source label"] = "Directorate National Security"
df["Notes"] = pd.NA

df['Date'] = pd.to_datetime(df['Date'].astype(str)).dt.date

df.to_csv("Suriname.csv", index=False)
