import csv
import json
from urllib.request import urlopen

import pandas as pd


# Data Sources
# https://outbreak.info/situation-reports/methods
# https://cov-lineages.org/index.html
# https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/


def get_locations():
    response = urlopen("https://api.outbreak.info/genomics/location?name=**")
    json_data = response.read().decode('utf-8', 'replace')
    loc_json = json.loads(json_data)

    loc_df = pd.json_normalize(loc_json["results"])
    loc_df = loc_df[["country_id", "country"]].drop_duplicates()
    loc_df = loc_df[loc_df["country_id"] != "None"]
    loc_df = loc_df.sort_values(by=['country'])

    return loc_df


def get_location_data(location_id, location):
    json_data = urlopen(
        f"https://api.outbreak.info/genomics/prevalence-by-location-all-lineages?location_id={location_id}&"
        f"cumulative=false&other_threshold=0.0&nday_threshold=0&ndays=1024").read().decode('utf-8', 'replace')
    loc_df = pd.json_normalize(json.loads(json_data)["results"])
    loc_df["location"] = location
    return loc_df


def main():
    locations_list = []

    locations = get_locations().to_dict('records')

    for location in locations:
        print(f"Location: {location['country']}")
        df = get_location_data(location["country_id"], location["country"])
        df.rename(
            columns={'lineage': 'variant', 'prevalence_rolling': 'perc_sequences', 'total_count': 'num_sequences_total',
                     'lineage_count': 'num_sequences'}, inplace=True)
        df.drop(['prevalence'], axis=1)

        df = df.reindex(
            columns=['location', 'date', 'variant', 'num_sequences', 'perc_sequences', 'num_sequences_total'])

        locations_list.append(df)

    df = pd.concat(locations_list)

    lineage_map = {
        # VOC
        "^B\\.1\\.1\\.7(.*)": "Alpha (VOC)",
        "^B\\.1\\.351(.*)": "Beta (VOC)",
        "^P\\.1(.*)": "Gamma - P.1 (VOC)",
        "^B\\.1\\.617\\.2": "Delta (VOC)",
        # VOI
        "^B\\.1\\.617\\.1": "Kappa (VOI)",
        "^B\\.1\\.525": "Eta (VOI)",
        "^B\\.1\\.526": "Iota (VOI)",
        "^C\\.37": "Lambda (VOI)",
        # Alerts for Further Monitoring
        "^B\\.1\\.427": "Epsilon (AFM)",
        "^B\\.1\\.429": "Epsilon (AFM)",
        "^P\\.2(.*)": "Zeta - P.2 (AFM)",
        "^P\\.3(.*)": "Theta - P.3 (AFM)",
        "^R\\.1(.*)": "R.1 (AFM)",
        "^R\\.2(.*)": "R.2 (AFM)",
        "^B\\.1\\.466\\.2": "B.1.466.2 (AFM)",
        "^B\\.1\\.621": "B.1.621 (AFM)",
        "^AV\\.1": "AV.1 (AFM)",
        "^B\\.1\\.1\\.318": "B.1.1.318 (AFM)",
        "^B\\.1\\.1\\.519": "B.1.1.519 (AFM)",
        "^AT\\.1": "AT.1 (AFM)",
        "^C\\.36\\.3": "C.36.3 (AFM)",
        "^C\\.36\\.3\\.1": "C.36.3.1 (AFM)",
        "^B\\.1\\.214\\.2": "B.1.214.2 (AFM)",
        # Others
        "^B\\.1\\.177(.*)": "B.1.177 (20E/EU1)",
        "^B\\.1\\.1\\.28": "B.1.1.28",
        "^OTHER": "Other"
    }

    # Clear zeroes
    df = df[df.perc_sequences != 0]

    df['variant'] = df['variant'].str.upper()

    df["variant"].replace(lineage_map, inplace=True, regex=True)

    main_lineage = [v for k, v in lineage_map.items()]
    other_lineage = [l for l in df["variant"].unique() if l not in main_lineage]
    lineage_to_parent = {}
    for o in other_lineage:
        lineage_to_parent[o] = o.split(".")[0] + " (Lineage)"

    df["variant"].replace(lineage_to_parent, inplace=True, regex=False)

    df['date'] = pd.to_datetime(df['date'])
    df = df.groupby(['location', pd.Grouper(key='date', freq='2W'), 'variant']).mean().reset_index().sort_values(
        ['location', 'date'])

    df.to_csv("../data/genomics.csv", index=False, quoting=csv.QUOTE_ALL, decimal=",")


main()
