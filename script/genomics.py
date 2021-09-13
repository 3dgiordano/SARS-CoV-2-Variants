import csv
import json
import re
from datetime import datetime, timedelta
from urllib.request import urlopen

import pandas as pd

# Data Sources
# https://outbreak.info/situation-reports/methods
# https://cov-lineages.org/index.html
# https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/
# https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html
# https://www.ecdc.europa.eu/en/covid-19/variants-concern
# https://www.gov.uk/government/collections/new-sars-cov-2-variant

interest_map = {"WHO": 1000, "CDC": 2000, "ECDC": 3000, "UY-GTI": 4000, "": 9000}
interest_type_map = {"VOC": 100, "VOI": 200, "AFM": 300, "VUM": 400, "": 900}

lineage_map = {
    "^A\\.23\\.1(.*)": "A.23.1",
    "^A\\.27(.*)": "A.27",
    "^A\\.28(.*)": "A.28",
    "^B\\.1\\.160(.*)": "B.1.160 - 20A/EU2",
    "^B\\.1\\.177(.*)": "B.1.177 - 20E/EU1",
    "^B\\.1\\.221(.*)": "B.1.221 - 20A/S:98F",
    "^B\\.1\\.258(.*)": "B.1.258 - 20A/S:439K",
    "^B\\.1\\.367(.*)": "B.1.367 - 20C/S:80Y",
    "^B\\.1\\.620(.*)": "B.1.620 - 20A/S:126A",
    "^B\\.1\\.616(.*)": "B.1.616",
    "^B\\.1\\.671\\.2$": "B.1.671.2",
    "^B\\.1\\.1\\.28$": "B.1.1.28",
    "^B\\.1\\.1\\.277(.*)": "B.1.1.277 - 20B/S:626S",
    "^B\\.1\\.1\\.302(.*)": "B.1.1.302 - 20B/S:1122L",
    "^B\\.1\\.1\\.519(.*)": "B.1.1.519 - 20B/S:732A",
    "^C\\.16(.*)": "C.16",
    "^P\\.6(.*)": "P.6 (UY-GTI VOC)",
    "^P\\.7(.*)": "P.7",
    "^R\\.2(.*)": "R.2",
    "^Q\\.1(.*)": "Alpha - Q.1 (WHO VOC)",
    "^Q\\.2(.*)": "Alpha - Q.2 (WHO VOC)",
    "^Q\\.3(.*)": "Alpha - Q.3 (WHO VOC)",
    "^Q\\.4(.*)": "Alpha - Q.4 (WHO VOC)",
    "^Q\\.5(.*)": "Alpha - Q.5 (WHO VOC)",
    "^Q\\.6(.*)": "Alpha - Q.6 (WHO VOC)",
    "^Q\\.7(.*)": "Alpha - Q.7 (WHO VOC)",
    "^Q\\.8(.*)": "Alpha - Q.8 (WHO VOC)",

    "^AY\\.1(.*)": "Delta - AY.1 (WHO VOC)",
    "^AY\\.2(.*)": "Delta - AY.2 (WHO VOC)",
    "^AY\\.3$": "Delta - AY.3 (WHO VOC)",
    "^AY\\.3\\.1$": "Delta - AY.3.1 (WHO VOC)",
    "^AY\\.4(.*)": "Delta - AY.4 (WHO VOC)",
    "^AY\\.5(.*)": "Delta - AY.5 (WHO VOC)",
    "^AY\\.6(.*)": "Delta - AY.6 (WHO VOC)",
    "^AY\\.7(.*)": "Delta - AY.7 (WHO VOC)",
    "^AY\\.8(.*)": "Delta - AY.8 (WHO VOC)",
    "^AY\\.9(.*)": "Delta - AY.9 (WHO VOC)",
    "^AY\\.10(.*)": "Delta - AY.10 (WHO VOC)",
    "^AY\\.11(.*)": "Delta - AY.11 (WHO VOC)",
    "^AY\\.12(.*)": "Delta - AY.12 (WHO VOC)",
    "^AY\\.13(.*)": "Delta - AY.13 (WHO VOC)",
    "^AY\\.14(.*)": "Delta - AY.14 (WHO VOC)",
    "^AY\\.15(.*)": "Delta - AY.15 (WHO VOC)",
    "^AY\\.16(.*)": "Delta - AY.16 (WHO VOC)",
    "^AY\\.17(.*)": "Delta - AY.17 (WHO VOC)",
    "^AY\\.18(.*)": "Delta - AY.18 (WHO VOC)",
    "^AY\\.19(.*)": "Delta - AY.19 (WHO VOC)",
    "^AY\\.20(.*)": "Delta - AY.20 (WHO VOC)",
    "^AY\\.21(.*)": "Delta - AY.21 (WHO VOC)",
    "^AY\\.22(.*)": "Delta - AY.22 (WHO VOC)",
    "^AY\\.23(.*)": "Delta - AY.23 (WHO VOC)",
    "^AY\\.24(.*)": "Delta - AY.24 (WHO VOC)",
    "^AY\\.25(.*)": "Delta - AY.25 (WHO VOC)",
    "^AY\\.26(.*)": "Delta - AY.26 (WHO VOC)",
    "^AY\\.27(.*)": "Delta - AY.27 (WHO VOC)",
    "^AY\\.28(.*)": "Delta - AY.28 (WHO VOC)",
    "^AY\\.29(.*)": "Delta - AY.29 (WHO VOC)",
    "^AY\\.30(.*)": "Delta - AY.30 (WHO VOC)",
    "^AY\\.31(.*)": "Delta - AY.31 (WHO VOC)",
    "^AY\\.32(.*)": "Delta - AY.32 (WHO VOC)",
    "^OTHER$": "Other"
}

who_detail_map = {
    "Gamma": "Gamma - P.1",
    "Mu": "Mu - 21H",
}

who_pango_map = {
    "^B\\.1\\.427(.*)": "Epsilon",
    "^B\\.1\\.429(.*)": "Epsilon",
    "^P\\.2(.*)": "Zeta - P.2",
    "^P\\.3(.*)": "Theta - P.3",
}


def get_url(url):
    from urllib.error import URLError
    tries = 3
    timeout = 60000
    for _ in range(tries):
        try:
            response = urlopen(url, timeout=timeout)
            break
        except URLError as err:
            if not isinstance(err.reason, socket.timeout):
                raise
    else:
        raise err
    return response


def get_locations():
    response = get_url("https://api.outbreak.info/genomics/location?name=**")
    json_data = response.read().decode('utf-8', 'replace')
    loc_json = json.loads(json_data)

    loc_df = pd.json_normalize(loc_json["results"])
    loc_df = loc_df[["country_id", "country"]].drop_duplicates()
    loc_df = loc_df[loc_df["country_id"] != "None"]
    loc_df = loc_df.sort_values(by=['country'])

    return loc_df


def get_location_data(location_id, location):
    json_data = get_url(
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
    return data.assign(pango=data['Pango lineage*'].str.split()).explode('pango')


def get_who_variants():
    from urllib.request import Request, urlopen

    who_variants_tracking_url = "https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/"

    who_body = urlopen(Request(who_variants_tracking_url, headers={'User-Agent': 'Mozilla/5.0'})).read().decode('UTF-8')

    # WA: generate spaces to solve new lines and unexpected new lines with div/p inside column values
    who_body = who_body.replace("<br />", "<br />&nbsp;").replace("</div>", "</div>&nbsp;").replace("</p>",
                                                                                                    "</p>&nbsp;")
    # Workaround
    who_body = who_body.replace("#", "")
    who_body = who_body.replace('<sup>&sect;</sup>', "")
    # who_body = who_body.replace("B.1.1.7", "B.1.1.7 Q ")
    # who_body = who_body.replace("B.1.617.2", "B.1.617.2 AY ")

    return pd.read_html(who_body, match=r'GISAID\sclade')


def cdc_filter_variants(table):
    pango_list = []
    for row in table.find_all(role='row'):
        for el in row.find_all(role="cell")[0].find_all('p'):
            text = el.text
            if text.startswith("Pango"):
                pango_list += [v.strip() for v in
                               text.split(":")[1].split("(")[0].strip().split(" and ")[0].split(",")]
    return pango_list


def get_cdc_variants():
    from urllib.request import Request, urlopen
    from bs4 import BeautifulSoup

    cdc_variants_tracking_url = "https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html"

    cdc_body = urlopen(Request(cdc_variants_tracking_url, headers={'User-Agent': 'Mozilla/5.0'})).read().decode('UTF-8')

    soup = BeautifulSoup(cdc_body, 'lxml')
    variants_tables = soup.find_all(
        'div', role="table")

    return cdc_filter_variants(variants_tables[0]), cdc_filter_variants(variants_tables[1])


def filter_to_dict(who_dict_map, cdc_variants, cdc_type):
    variants = {}
    for cdc_variant in cdc_variants:
        match = False
        for k, v in who_dict_map.items():
            if re.match(k, cdc_variant):
                match = True
                break
        if not match:
            variant_regex = pango_regex(cdc_variant)
            # Invert, check if this match on who variant, if match try to fix on that lineage
            match = False
            k_normal = ""
            for k, v in who_dict_map.items():
                k_normal = k.replace("\\", "").replace("^", "").replace("$", "").replace("(.*)", "")
                if re.match(variant_regex, k_normal):
                    match = True
                    break
            if match:
                variant_regex = variant_regex.replace("(.*)", "$")
            if k_normal != variant_regex:
                variants[variant_regex] = f"{who_pango_rename(cdc_variant)} {cdc_type}"

    return variants


def ecdc_filter_values(table):
    return list(set([x.split("+")[0].strip() for x in table['Lineage + additional mutations'].tolist()]))


def get_ecdc_variants():
    from urllib.request import Request, urlopen

    ecdc_variants_tracking_url = "https://www.ecdc.europa.eu/en/covid-19/variants-concern"

    ecdc_body = urlopen(Request(ecdc_variants_tracking_url, headers={'User-Agent': 'Mozilla/5.0'})).read().decode(
        'UTF-8')

    ecdc_tables = pd.read_html(ecdc_body, match=r'WHO\slabel')

    ecdc_voc = ecdc_filter_values(ecdc_tables[0])
    ecdc_voi = list(set(ecdc_filter_values(ecdc_tables[1])) - set(ecdc_voc))
    ecdc_vum = list((set(ecdc_filter_values(ecdc_tables[2])) - set(ecdc_voi)) - set(ecdc_voi))
    return ecdc_voc, ecdc_voi, ecdc_vum


def phe_filter_values(table, v_type):
    # WA wrong pango from B.1.324.1 to B.1.623
    return list(set(
        [x.split("PANGO: ")[1].split("nextstrain: ")[0].replace("B.1.324.1", "B.1.623").split(",")[0].strip() for x in
         table.loc[table["Label"].str.startswith(v_type, na=False)]['Lineages'].tolist() if
         "Multiple" not in x]))


def get_phe_variants():
    from urllib.request import Request, urlopen

    phe_variants_tracking_url = "https://github.com/phe-genomics/variant_definitions/blob/main/README.md"

    phe_body = urlopen(Request(phe_variants_tracking_url, headers={'User-Agent': 'Mozilla/5.0'})).read().decode('UTF-8')

    phe_tables = pd.read_html(phe_body, match=r'Lineages')

    phe_voc = phe_filter_values(phe_tables[0], 'VOC')
    phe_vui = phe_filter_values(phe_tables[0], 'VUI')

    return phe_voc, phe_vui


def get_lineage_map():
    (who_voc, who_voi, who_afm) = get_who_variants()
    (cdc_voi, cdc_voc) = get_cdc_variants()
    (ecdc_voc, ecdc_voi, ecdc_vum) = get_ecdc_variants()
    (phe_voc, phe_vui) = get_phe_variants()

    who_voc = who_expand(who_voc)
    who_voi = who_expand(who_voi)
    who_afm = who_expand(who_afm)

    who_voc_dict = who_to_dict(who_voc, "(WHO VOC)")
    who_voi_dict = who_to_dict(who_voi, "(WHO VOI)")
    who_afm_dict = who_to_dict(who_afm, "(WHO AFM)")

    lineage_dict_map = dict(
        list(who_voc_dict.items()) +
        list(who_voi_dict.items()) +
        list(who_afm_dict.items())
    )

    cdc_voc_dict = filter_to_dict(lineage_dict_map, cdc_voc, "(CDC VOC)")
    cdc_voi_dict = filter_to_dict(lineage_dict_map, cdc_voi, "(CDC VOI)")

    lineage_dict_map.update(cdc_voi_dict)
    lineage_dict_map.update(cdc_voc_dict)

    ecdc_voc_dict = filter_to_dict(lineage_dict_map, ecdc_voc, "(ECDC VOC)")
    ecdc_voi_dict = filter_to_dict(lineage_dict_map, ecdc_voi, "(ECDC VOI)")
    ecdc_vum_dict = filter_to_dict(lineage_dict_map, ecdc_vum, "(ECDC VUM)")

    lineage_dict_map.update(ecdc_voi_dict)
    lineage_dict_map.update(ecdc_voc_dict)
    lineage_dict_map.update(ecdc_vum_dict)

    phe_voc_dict = filter_to_dict(lineage_dict_map, phe_voc, "(PHE VOC)")
    phe_vui_dict = filter_to_dict(lineage_dict_map, phe_vui, "(PHE VUI)")

    lineage_dict_map.update(phe_voc_dict)
    lineage_dict_map.update(phe_vui_dict)

    lineage_dict_map.update(lineage_map)

    data = {
        "map": lineage_dict_map,
        "data": {
            "who": {"voc": who_voc, "voi": who_voc, "afm": who_afm},
            "cdc": {"voi": cdc_voi, "voc": cdc_voc},
            "ecdc": {"voi": ecdc_voi, "voc": ecdc_voc, "vum": ecdc_vum},
            "phe": {"voc": phe_voc, "vui": phe_vui},
        }
    }
    return data


def get_sort_order(interest, interest_type, label):
    if interest not in interest_map:
        interest = ""
    if interest_type not in interest_type_map:
        interest_type = ""
    return str(interest_map[interest] + interest_type_map[interest_type]) + label


def export_variants(main_lineage_map):
    main_lineage = list(dict.fromkeys([v for k, v in main_lineage_map.items()]))
    variant_records = []

    for lineage in main_lineage:
        if lineage == "Other":
            continue

        label = lineage.split("(")[0].strip()
        pango_list = ", ".join([k.replace("\\", "").replace("(.*)", ".*").replace("^", "").replace("$", "") for k, v in
                                main_lineage_map.items() if v == lineage])

        interest = lineage.split("(")[1].replace(")", "").split()[0] if \
            lineage.find("(") > -1 else "" if ["VUM", "VOI", "VOC", "AFM"] else ""

        interest_type = ""
        if "(" in lineage and len(lineage.split("(")[1].replace(")", "").split()) > 0:
            interest_type = lineage.split("(")[1].replace(")", "").split()[1]

        variant_records.append({
            "label": label,
            "pango": pango_list,
            "interest": interest,
            "type": interest_type,
            "sort_order": get_sort_order(interest, interest_type, label)
        })
    pd.json_normalize(variant_records).sort_values(["sort_order"]).drop(columns=['sort_order']).to_csv(
        f"../data/variants.csv", index=False,
        quoting=csv.QUOTE_ALL, decimal=",")


def main():
    locations_list = []

    locations = get_locations().to_dict('records')

    for location in locations:

        # if location["country"] != "Uruguay":
        #     continue

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

    data = get_lineage_map()
    main_lineage_map = data["map"]

    df["variant"].replace(main_lineage_map, inplace=True, regex=True)

    main_lineage = list(dict.fromkeys([v for k, v in main_lineage_map.items()]))
    other_lineage = list(dict.fromkeys([l for l in df["variant"].unique() if l not in main_lineage]))

    export_variants(main_lineage_map)

    # Transform no specific lineage to parent lineage
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
    # Exclude the last register because is noisy
    to_date = datetime.now() - timedelta(days=14)
    for location in locations:
        print(f"Save Location: {location['country']}")
        df_location = df_pivoted[df_pivoted["location"] == location['country']]
        df_location = df_location.loc[:, (df_location != 0).any(axis=0)]  # Remove zeroes columns
        if df_location.size > 0:
            if df_location.tail(1).iloc[0]['date'] > to_date:
                # Remove the last register because probably is noise
                df_location = df_location[:-1]
        df_location.to_csv(
            f"../data/{location['country']}.csv", index=False, quoting=csv.QUOTE_ALL, decimal=",")

    df_world_pivoted = df.groupby(['date', 'variant']).agg(
        {'perc_sequences': 'mean'}).reset_index().pivot(index=["date"], columns=["variant"],
                                                        values="perc_sequences").reset_index()

    df_world_pivoted["location"] = "World"
    df_world_pivoted.insert(0, "location", df_world_pivoted.pop("location"))
    df_world_pivoted = df_world_pivoted.fillna(0)

    df_world_pivoted[:-3].to_csv("../data/World.csv", index=False, quoting=csv.QUOTE_ALL, decimal=",")

    update_dict = {
        "last_update_utc": datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S.%f%z')
    }
    with open('../data/update.json', 'w') as upd_file:
        json.dump(update_dict, upd_file, indent=4)


main()
