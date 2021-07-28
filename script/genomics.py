import csv
import json
import re
from urllib.request import urlopen

import pandas as pd

# Data Sources
# https://outbreak.info/situation-reports/methods
# https://cov-lineages.org/index.html
# https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/

lineage_map = {
    "^B\\.1\\.177(.*)": "B.1.177 (20E/EU1)",
    "^B\\.1\\.1\\.28$": "B.1.1.28",
    "^OTHER$": "Other"
}

who_detail_map = {
    "Gamma": "Gamma - P.1",
}

who_pango_map = {
    "^B\\.1\\.427(.*)": "Epsilon",
    "^B\\.1\\.429(.*)": "Epsilon",
    "^P\\.2(.*)": "Zeta - P.2",
    "^P\\.3(.*)": "Theta - P.3"
}


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


def who_detail(who_name):
    for k, v in who_detail_map.items():
        who_name = who_name.replace(k, v)
    return who_name


def who_pango_rename(pango):
    for k, v in who_pango_map.items():
        pango = re.sub(k, v, pango)
    return pango


def pango_regex(pango):
    sub_depth = pango.count(".")
    pango = pango.replace("*", "").replace(".", "\\.")
    return f"^{pango}$" if sub_depth == 3 else f"^{pango}(.*)"


def who_to_dict(data, who_type):
    who_label = 'WHO\xa0label'
    if who_label in data.columns:
        return {pango_regex(row['pango']): f"{who_detail(row[who_label])} {who_type}" for
                ind, row in
                data.iterrows()}
    else:
        return {
            pango_regex(row['pango']):
                f"{who_pango_rename(row['pango'].replace('*', ''))} {who_type}" for ind, row in
            data.iterrows()}


def who_expand(data):
    return data.assign(pango=data['Pango lineages'].str.split()).explode('pango')


def get_lineage_map():
    from urllib.request import Request, urlopen

    who_variants_tracking_url = "https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/"

    who_body = urlopen(Request(who_variants_tracking_url, headers={'User-Agent': 'Mozilla/5.0'})).read().decode('UTF-8')

    # WA: generate spaces to solve new lines and unexpected new lines with div/p inside column values
    who_body = who_body.replace("<br />", "<br />&nbsp;").replace("</div>", "</div>&nbsp;").replace("</p>",
                                                                                                    "</p>&nbsp;")

    (who_voc, who_voi, who_afm) = pd.read_html(who_body, match=r'GISAID\sclade')

    who_voc = who_expand(who_voc)
    who_voi = who_expand(who_voi)
    who_afm = who_expand(who_afm)

    who_voc_dict = who_to_dict(who_voc, "(VOC)")
    who_voi_dict = who_to_dict(who_voi, "(VOI)")
    who_afm_dict = who_to_dict(who_afm, "(AFM)")

    who_dict_map = dict(
        list(who_voc_dict.items()) +
        list(who_voi_dict.items()) +
        list(who_afm_dict.items())
    )
    who_dict_map.update(lineage_map)

    return who_dict_map


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

    # Clear zeroes
    df = df[df.perc_sequences != 0]

    df['variant'] = df['variant'].str.upper()

    who_lineage_map = get_lineage_map()

    df["variant"].replace(who_lineage_map, inplace=True, regex=True)

    main_lineage = list(dict.fromkeys([v for k, v in who_lineage_map.items()]))
    other_lineage = list(dict.fromkeys([l for l in df["variant"].unique() if l not in main_lineage]))

    lineage_to_parent = {}
    for o in other_lineage:
        lineage_to_parent[o] = o.split(".")[0] + " (Lineage)"

    df["variant"].replace(lineage_to_parent, inplace=True, regex=False)

    df['date'] = pd.to_datetime(df['date'])

    df = df.groupby(['location', pd.Grouper(key='date', freq='2W'), 'variant']).agg(
        {'num_sequences': 'sum'}).reset_index()

    dfb = df.groupby(['location', 'date']).agg(
        {'num_sequences': 'sum'}).rename(columns={"num_sequences": "num_sequences_total"}).reset_index()

    df = pd.merge(df, dfb, on=['location', 'date'])

    df["perc_sequences"] = (df["num_sequences"] / df["num_sequences_total"]) * 100

    df = df.drop(df[df.perc_sequences.isnull()].index)

    df = df.sort_values(['location', 'date', 'variant'])

    df.to_csv("../data/genomics.csv", index=False, quoting=csv.QUOTE_ALL, decimal=",")

    df_pivoted = df.pivot(index=["location", "date"], columns=["variant"], values="perc_sequences").reset_index()
    df_pivoted = df_pivoted.fillna(0)

    # Save a file for each location generating pivot table
    for location in locations:
        print(f"Save Location: {location['country']}")
        df_pivoted[df_pivoted["location"] == location['country']].to_csv(f"../data/{location['country']}.csv",
                                                                         index=False,
                                                                         quoting=csv.QUOTE_ALL, decimal=",")


main()
