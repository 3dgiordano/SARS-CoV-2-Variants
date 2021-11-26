import csv
import json
import os
import re
from datetime import datetime, timedelta
from urllib.request import urlopen

import math
import numpy as np
import pandas as pd
from numpy import errstate, isneginf, array

# Data Sources
# https://outbreak.info/situation-reports/methods
# https://cov-lineages.org/index.html
# https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/
# https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html
# https://www.ecdc.europa.eu/en/covid-19/variants-concern
# https://www.gov.uk/government/collections/new-sars-cov-2-variant

interest_map = {"WHO": 1000, "CDC": 2000, "ECDC": 3000, "UY-GTI": 4000, "": 9000}
interest_type_map = {"VOC": 100, "VOI": 200, "AFM": 300, "VUM": 400, "": 900}
lineages = None
t = "($|\\..*$)"  # Tail Regex
lineage_map = {
    f"^A\\.2\\.5{t}": "A.2.5",
    f"^A\\.23\\.1{t}": "A.23.1",
    f"^A\\.27{t}": "A.27",
    f"^A\\.28{t}": "A.28",
    f"^AT\\.1{t}": "AT.1",
    f"^B\\.1\\.160{t}": "B.1.160 - 20A/EU2",
    f"^B\\.1\\.177{t}": "B.1.177 - 20E/EU1",
    f"^B\\.1\\.221{t}": "B.1.221 - 20A/S:98F",
    f"^B\\.1\\.258{t}": "B.1.258 - 20A/S:439K",
    f"^B\\.1\\.367{t}": "B.1.367 - 20C/S:80Y",
    f"^B\\.1\\.620{t}": "B.1.620 - 20A/S:126A",
    f"^B\\.1\\.616{t}": "B.1.616",
    f"^B\\.1\\.628{t}": "B.1.628",
    f"^B\\.1\\.671\\.2$": "B.1.671.2",
    f"^B\\.1\\.1\\.28$": "B.1.1.28",
    f"^B\\.1\\.1\\.277{t}": "B.1.1.277 - 20B/S:626S",
    f"^B\\.1\\.1\\.302{t}": "B.1.1.302 - 20B/S:1122L",
    f"^B\\.1\\.1\\.519{t}": "B.1.1.519 - 20B/S:732A",
    f"^C\\.16{t}": "C.16",
    f"^N\\.7{t}": "N.7 (UY-GTI)",
    f"^P\\.6{t}": "P.6 (UY-GTI)",
    f"^P\\.7{t}": "P.7",
    f"^R\\.2{t}": "R.2",
    "^OTHER$": "Other"
}

who_detail_map = {
    # "Gamma": "Gamma - P.1",
    "Mu": "Mu - 21H",
}

who_pango_map = {
    f"^B\\.1\\.427{t}": "Epsilon",
    f"^B\\.1\\.429{t}": "Epsilon",
    f"^P\\.1{t}": "Gamma - P.1",
    f"^P\\.2{t}": "Zeta - P.2",
    f"^P\\.3{t}": "Theta - P.3",
    f"^B\\.1\\.525{t}": "Eta",
    f"^B\\.1\\.526{t}": "Iota",
    f"^B\\.1\\.617\\.1$": "Kappa",
}


def get_url(url):
    from urllib.error import URLError
    import socket
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


def get_cases_data():
    return pd.read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/jhu/full_data.csv")


def get_locations_data():
    return pd.read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/jhu/locations.csv")


def get_cases_r_data():
    r_file = f"../temp/r_data.csv"
    last_time = None
    if os.path.isfile(r_file):
        last_time = (datetime.now() - datetime.fromtimestamp(os.path.getmtime(r_file))).total_seconds()
    # keep the stored data temporarily for one hour
    if not last_time or last_time > 21600:
        r_df = pd.read_csv("https://raw.githubusercontent.com/crondonm/TrackingR/main/Estimates-Database/database.csv")
        r_df.to_csv(r_file, index=False, quoting=csv.QUOTE_ALL, decimal=",")
    else:
        r_df = pd.read_csv(r_file, quoting=csv.QUOTE_ALL, decimal=",")
    return r_df


def get_alias_map_sub_lineage(lineage_to_match):
    global lineages
    alias_map_sub_lineage = []
    if not lineages:
        response = get_url("https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt")
        lineages = response.read().decode('utf-8', 'replace').splitlines()
    for lineage in lineages:
        alias = lineage.split("\t")[0]
        if "\t" in lineage:
            is_alias_of = lineage.split("\t")[1].split("Alias of ")
        else:
            is_alias_of = lineage.split(", ")[1].split("Alias of ")
        if len(is_alias_of) > 1 and not alias.startswith("*"):
            alias_of = is_alias_of[1].split(",")[0]
            alias_of_x = ".".join(alias_of.split(".")[:lineage_to_match.count(".") + 1])
            if lineage_to_match == alias_of_x:
                alias_map_sub_lineage.append(alias)
    return alias_map_sub_lineage


def get_locations():
    locations_file = f"../temp/locations_data.csv"
    last_time = None
    if os.path.isfile(locations_file):
        last_time = (datetime.now() - datetime.fromtimestamp(os.path.getmtime(locations_file))).total_seconds()
    # keep the stored data temporarily for one hour
    if not last_time or last_time > 21600:
        response = get_url("https://api.outbreak.info/genomics/location?name=**")
        json_data = response.read().decode('utf-8', 'replace')
        loc_json = json.loads(json_data)

        loc_df = pd.json_normalize(loc_json["results"])
        loc_df = loc_df[["country_id", "country"]].drop_duplicates()
        loc_df = loc_df[loc_df["country_id"] != "None"]
        loc_df = loc_df.sort_values(by=['country'])
        loc_df.to_csv(locations_file, index=False, quoting=csv.QUOTE_ALL, decimal=",")
    else:
        loc_df = pd.read_csv(locations_file, quoting=csv.QUOTE_ALL, decimal=",")

    return loc_df


def get_location_data(location_id, location):
    location_file = f"../temp/{location}_data.csv"
    last_time = None
    if os.path.isfile(location_file):
        last_time = (datetime.now() - datetime.fromtimestamp(os.path.getmtime(location_file))).total_seconds()
    # keep the stored data temporarily for one hour
    if not last_time or last_time > 21600:
        json_data = get_url(
            f"https://api.outbreak.info/genomics/prevalence-by-location-all-lineages?location_id={location_id}&"
            f"cumulative=false&other_threshold=0.0&nday_threshold=0&ndays=2048").read().decode('utf-8', 'replace')
        loc_df = pd.json_normalize(json.loads(json_data)["results"])
        loc_df["location"] = location
        loc_df = loc_df.sort_values(['location', 'date', 'lineage'])
        loc_df.to_csv(location_file, index=False, quoting=csv.QUOTE_ALL, decimal=",")
    else:
        loc_df = pd.read_csv(location_file, quoting=csv.QUOTE_ALL, decimal=",")
        # loc_df.convert_objects(convert_numeric=True)
    return loc_df


def who_detail(who_name):
    for k, v in who_detail_map.items():
        who_name = who_name.replace(k, v)
    return who_name


def who_pango_rename(pango):
    for k, v in who_pango_map.items():
        pango = re.sub(k, v, pango)

    return pango


def pango_regex(pango, no_sub=False):
    sub_depth = pango.count(".")
    pango = pango.replace("*", "").replace(".", "\\.")
    return f"^{pango}$" if sub_depth == 3 or no_sub else f"^{pango}{t}"


def who_to_dict(data, who_type):
    who_label = 'WHO\xa0label'
    with_who_label = False
    if who_label in data.columns:
        with_who_label = True

    who_dict = {}
    for ind, row in data.iterrows():
        pango = row['pango']

        if with_who_label:
            label = who_detail(row[who_label])
        else:
            label = f"{who_pango_rename(row['pango'].replace('*', ''))}"

        if pango == "P.1":
            who_dict["^P.1$"] = f"{label} {who_type}"
            pango_alias_lineages = get_alias_map_sub_lineage("B.1.1.28.1")
        else:
            who_dict[pango_regex(pango)] = f"{label} {who_type}"
            pango_alias_lineages = get_alias_map_sub_lineage(pango)
        if pango_alias_lineages:
            for p_alias_l in pango_alias_lineages:
                who_dict[pango_regex(p_alias_l, no_sub=True)] = f"{label} - {p_alias_l} {who_type}"

    return who_dict


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
    who_body = who_body.replace("Omicron*", "Omicron")
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


def cdc_to_dict(data, who_type):
    who_label = 'WHO\xa0Label'
    if who_label in data.columns:
        with_who_label = True
    else:
        with_who_label = False
    who_dict = {}
    for ind, row in data.iterrows():
        pango = row['pango']
        if with_who_label:
            label = who_detail(row[who_label])
        else:
            label = f"{who_pango_rename(row['pango'].replace('*', ''))}"
        who_dict[pango_regex(pango)] = f"{label} {who_type}"
        pango_alias_lineages = get_alias_map_sub_lineage(pango)
        if pango_alias_lineages:
            for p_alias_l in pango_alias_lineages:
                who_dict[pango_regex(p_alias_l, no_sub=True)] = f"{label} - {p_alias_l} {who_type}"

    return who_dict


def cdc_expand(data):
    df = data.assign(pangoy=data['Pango Lineage'].str.split(",")).explode('pangoy')
    df["pangoy"] = df["pangoy"].str.strip()
    df = df.assign(pango=df["pangoy"].str.split()).explode("pango")
    return df


def get_cdc_variants():
    from urllib.request import Request, urlopen
    from bs4 import BeautifulSoup

    cdc_variants_tracking_url = "https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html"

    cdc_body = urlopen(Request(cdc_variants_tracking_url, headers={'User-Agent': 'Mozilla/5.0'})).read().decode('UTF-8')

    cdc_body = cdc_body.replace("and descendent lineages", "").replace(" and Q lineages", "")

    soup = BeautifulSoup(cdc_body, 'lxml')
    variants_tables = soup.find_all(
        'div', role="table")

    cdc_data = pd.read_html(cdc_body, match=r'Pango\sLineage')

    cdc_voi = []  # TODO:For now, CDC not have VOI variants
    return cdc_voi, cdc_filter_variants(variants_tables[0]), cdc_expand(cdc_data[0])["pango"].tolist()


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
                k_normal = k.replace("\\", "").replace("^", "").replace("$", "").replace(f"{t}", "")
                if re.match(variant_regex, k_normal):
                    match = True
                    break
            if match:
                variant_regex = variant_regex.replace(f"{t}", "$")
            if k_normal != variant_regex:
                variants[variant_regex] = f"{who_pango_rename(cdc_variant)} {cdc_type}"

    return variants


def ecdc_filter_values(table):
    return list(set([str(x).split("+")[0].split("(")[0].strip() for x in table['Lineage + additional mutations'].tolist()]))


def get_ecdc_variants():
    from urllib.request import Request, urlopen

    ecdc_variants_tracking_url = "https://www.ecdc.europa.eu/en/covid-19/variants-concern"

    ecdc_body = urlopen(Request(ecdc_variants_tracking_url, headers={'User-Agent': 'Mozilla/5.0'})).read().decode(
        'UTF-8')

    # Fix errors in HTML
    ecdc_body = re.sub(r"<td>\s+<table>\s+<tbody>\s+<tr>|</tr>\s+</tbody>\s+</table>\s+</td>", "", ecdc_body)

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
    (who_voc, who_voi, who_afm, who_fmv) = get_who_variants()
    (cdc_voi, cdc_voc, cdc_vbm) = get_cdc_variants()
    (ecdc_voc, ecdc_voi, ecdc_vum) = get_ecdc_variants()
    (phe_voc, phe_vui) = get_phe_variants()

    who_voc = who_expand(who_voc)
    who_voi = who_expand(who_voi)
    who_afm = who_expand(who_afm)
    who_fmv = who_expand(who_fmv)

    who_voc_dict = who_to_dict(who_voc, "(WHO VOC)")
    who_voi_dict = who_to_dict(who_voi, "(WHO VOI)")
    who_afm_dict = who_to_dict(who_afm, "(WHO AFM)")
    who_fmv_dict = who_to_dict(who_fmv, "(WHO FMV)")

    lineage_dict_map = dict(
        list(who_voc_dict.items()) +
        list(who_voi_dict.items()) +
        list(who_afm_dict.items()) +
        list(who_fmv_dict.items())
    )

    cdc_voc_dict = filter_to_dict(lineage_dict_map, cdc_voc, "(CDC VOC)")
    cdc_voi_dict = filter_to_dict(lineage_dict_map, cdc_voi, "(CDC VOI)")
    cdc_vbm_dict = filter_to_dict(lineage_dict_map, cdc_vbm, "(CDC VBM)")

    lineage_dict_map.update(cdc_voi_dict)
    lineage_dict_map.update(cdc_voc_dict)
    lineage_dict_map.update(cdc_vbm_dict)

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
            "who": {"voc": who_voc, "voi": who_voc, "afm": who_afm, "fmv": who_fmv},
            "cdc": {"voi": cdc_voi, "voc": cdc_voc, "vbm": cdc_vbm},
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
        pango_list = ", ".join([k.replace(f"{t}", ".*").replace("\\", "").replace("^", "").replace("$", "") for k, v in
                                main_lineage_map.items() if v == lineage])

        interest = lineage.split("(")[1].replace(")", "").split()[0] if \
            lineage.find("(") > -1 else "" if ["VUM", "VOI", "VOC", "AFM"] else ""

        interest_type = ""
        if "(" in lineage and len(lineage.split("(")[1].replace(")", "").split()) > 1:
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

    df_cases_data = get_cases_data()

    df_cases_data['date'] = pd.to_datetime(df_cases_data['date'])

    df_cases_data = df_cases_data.groupby(['location', pd.Grouper(key='date', freq='2W')]).agg(
        {'new_cases': 'sum'}).rename(columns={"new_cases": "cases"}).reset_index()

    df_cases_data["cases"] = pd.to_numeric(df_cases_data["cases"], downcast='integer')

    df_cases_data.cases = df_cases_data.cases.mask(df_cases_data.cases.lt(0), 0)  # Remove negative values

    for location in locations:
        # if location["country"] != "Uruguay":
        #     continue

        print(f"Location: {location['country']}")
        df = get_location_data(location["country_id"], location["country"])
        df.rename(
            columns={'lineage': 'variant', 'prevalence_rolling': 'perc_sequences', 'total_count': 'num_sequences_total',
                     'lineage_count': 'num_sequences'}, inplace=True)
        df.drop(['prevalence'], axis=1, inplace=True)

        df = df.reindex(
            columns=['location', 'date', 'variant', 'num_sequences', 'perc_sequences', 'num_sequences_total'])

        locations_list.append(df)

    print("Create location list...")
    df = pd.concat(locations_list)

    # Clear zeroes
    df = df[df.perc_sequences != 0]

    df['variant'] = df['variant'].str.upper()

    # Drop from cases countries not in list
    locations_to_clean = list(set(df_cases_data["location"].values) - set(df["location"].values))
    df_cases_data = df_cases_data.loc[~df_cases_data.location.isin(locations_to_clean)]

    print("Save cases.csv...")
    df_cases_data.to_csv("../data/cases.csv", index=False, quoting=csv.QUOTE_ALL, decimal=",")

    # R

    print("get looations")
    df_cases_location_data = get_locations_data()

    # Clean to only needed data
    df_cases_location_data = df_cases_location_data[
        df_cases_location_data.columns.intersection(["continent", "location", "population"])]

    print("get cases r")
    df_cases_r_data = get_cases_r_data()

    print("Rename countries")
    df_cases_r_data.rename(columns={"Country/Region": "location", "Date": "date", "R": "r"}, inplace=True)

    df_cases_r_data['date'] = pd.to_datetime(df_cases_r_data['date'])

    df_cases_r_data = df_cases_r_data.groupby(['location', pd.Grouper(key='date', freq='2W')]).agg(
        {'r': 'mean'}).reset_index()

    print("Merge cases")
    df_cases_r_data = pd.merge(df_cases_r_data, df_cases_data, on=['location', 'date'], how='outer', sort=True)
    df_cases_r_data["cases"] = df_cases_r_data["cases"].fillna(0)
    df_cases_r_data["cases"] = pd.to_numeric(df_cases_r_data["cases"], downcast='integer')

    # Fix border cases
    def fix_row_x(r):

        if pd.isnull(r["r"]):
            if r["location"] == df_cases_r_data.iloc[r.name - 1]["location"] and not pd.isnull(
                    df_cases_r_data.iloc[r.name - 1]["r"]):
                df_cases_r_data.loc[[r.name], "r"] = df_cases_r_data.iloc[r.name - 1]["r"]
            else:
                df_cases_r_data.loc[[r.name], "r"] = 0.0

    df_cases_r_data.apply(fix_row_x, axis=1)

    print("Merge location")
    df_cases_r_data = pd.merge(df_cases_r_data, df_cases_location_data, on=['location'])
    df_cases_r_data["population"] = df_cases_r_data["population"].fillna(0)
    df_cases_r_data["population"] = df_cases_r_data["population"].astype(int)

    def pro_cases(x):
        if x["cases"] > 0:
            days = (datetime.now() - x["date"]).days
            if days < 0:
                tot_days_data = (14 + days) - 1
                p_cases = round((x["cases"] * 14) / tot_days_data, 0)
                df_cases_r_data.loc[[x.name], "cases"] = p_cases

    df_cases_r_data.apply(pro_cases, axis=1)

    df_cases_r_data["cases_100k"] = round((df_cases_r_data["cases"] / df_cases_r_data["population"]) * 100000, 2)

    df_cases_r_data["risk"] = round((df_cases_r_data["cases_100k"] * (df_cases_r_data["r"] + 1)) / 250, 2)

    df_cases_r_data["x"] = 0  # init

    def row_x(x):
        prev_x = 0
        prev_r = 0
        prev_cases = 0
        inc_cases = 1
        if x.name > 0 and x["location"] == df_cases_r_data.iloc[x.name - 1]["location"]:
            prev_x = df_cases_r_data.iloc[x.name - 1]["x"]
            prev_r = df_cases_r_data.iloc[x.name - 1]["r"]
            prev_cases = df_cases_r_data.iloc[x.name - 1]["cases"]
            if prev_cases > 0:
                inc_cases = (x["cases"] / prev_cases)

        if x["r"] >= 0.85:
            if prev_r == 0:
                val_x = 0
            else:
                if x["r"] >= 1:
                    if x["r"] >= prev_r and prev_x > 0:
                        val_x = prev_x / 3
                    else:
                        val_x = 0.5
                else:
                    if prev_r < 0.85:
                        val_x = 0.25
                    elif x["r"] >= 0.95:
                        val_x = 0.15
                    else:
                        val_x = -(prev_x / 3)
            val = val_x + prev_x
            if val < 0:
                val = 0
            df_cases_r_data.loc[[x.name], "x"] = val
        else:
            if prev_r < 0.85:
                val_x = prev_x
            else:
                val_x = prev_x / 1.5
            val = prev_x - val_x

            if val < 0:
                val = 0
            df_cases_r_data.loc[[x.name], "x"] = val

        # days = (datetime.now() - x["date"]).days
        # if days < 0:
        #    tot_days_data = (14 + days) - 1
        #    pond = tot_days_data / 14
        #    if x["cases"] > 0:

        #        p_cases = (x["cases"] * 14) / tot_days_data
        #        diff = p_cases / x["cases"]
        #        dif_z = p_cases / prev_cases
        #        if dif_z > 1:
        #            dif_z = dif_z * 2
        #        else:
        #            dif_z = dif_z / 2
        #        xpp = ((df_cases_r_data.loc[[x.name], "x"].item() * (1 + pond)) * diff) * dif_z  # .item() * 14 ) / (14 + days)

        #        print(f"{str(x['location'])} X:{str(df_cases_r_data.loc[[x.name], 'x'].item())} C:{str(x['cases'])} PC:{str(p_cases)} TD:{tot_days_data} PO:{pond} DF:{diff} DZ:{dif_z} XPP:{str(xpp)}")

        #        # if prev_r >= 0.8 and x["r"] >= 1 and xpp > df_cases_r_data.loc[[x.name], "x"].item():
        #        #     df_cases_r_data.loc[[x.name], "x"] = xpp * x["r"]

        #        df_cases_r_data.loc[[x.name], "x"] = xpp

        #        df_cases_r_data.loc[[x.name], "cases_100k"] = round((p_cases / x["population"]) * 100000, 2)

    df_cases_r_data.apply(row_x, axis=1)

    with errstate(divide='ignore'):
        df_cases_r_data["riskX"] = np.exp(df_cases_r_data["cases_100k"] / 1000)
        df_cases_r_data.riskX = df_cases_r_data.riskX.mask(df_cases_r_data.riskX.lt(0.0), 0)  # Clean negatives values

    df_cases_r_data["risk2"] = round(
        (((df_cases_r_data["cases_100k"] * df_cases_r_data["riskX"]) * (1 + df_cases_r_data["x"])) / 50), 2)
    df_cases_r_data["risk3"] = np.sqrt(np.sqrt(df_cases_r_data["risk2"]))

    df_cases_r_data.risk3 = df_cases_r_data.risk3.mask(df_cases_r_data.risk3.lt(0.5), 0)  # Clean lowers values

    # df_cases_r_data["population"] = pd.to_numeric(df_cases_r_data["population"], downcast='integer')

    print("Save cases_r.csv...")
    df_cases_r_data.to_csv("../data/cases_r.csv", index=False, quoting=csv.QUOTE_ALL, decimal=",")

    # return

    print("Map lineage...")
    data = get_lineage_map()
    main_lineage_map = data["map"]

    df["variant"].replace(main_lineage_map, inplace=True, regex=True)

    main_lineage = list(dict.fromkeys([v for k, v in main_lineage_map.items()]))
    other_lineage = list(dict.fromkeys([str(l) for l in df["variant"].unique() if l not in main_lineage]))

    print("Save lineage map...")
    export_variants(main_lineage_map)

    print("Transform not monitored lineage to parent...")
    # Transform no specific lineage to parent lineage
    lineage_to_parent = {}
    for o in other_lineage:
        lineage_to_parent[o] = o.split(".")[0] + " (Lineage)"
    df["variant"].replace(lineage_to_parent, inplace=True, regex=False)

    print("Group and calculate percentage of sequences...")
    df['date'] = pd.to_datetime(df['date'])

    df = df.groupby(['location', pd.Grouper(key='date', freq='2W'), 'variant']).agg(
        {'num_sequences': 'sum'}).reset_index()

    dfb = df.groupby(['location', 'date']).agg(
        {'num_sequences': 'sum'}).rename(columns={"num_sequences": "num_sequences_total"}).reset_index()

    df = pd.merge(df, dfb, on=['location', 'date'])

    df["perc_sequences"] = (df["num_sequences"] / df["num_sequences_total"]) * 100

    df = df.drop(df[df.perc_sequences.isnull()].index)

    df = df.sort_values(['location', 'date', 'variant'])

    print("Save genomics.csv...")
    df.to_csv("../data/genomics.csv", index=False, quoting=csv.QUOTE_ALL, decimal=",")

    print("Pivot data...")
    df_pivoted = df.pivot(index=["location", "date"], columns=["variant"], values="perc_sequences").reset_index()
    df_pivoted = df_pivoted.fillna(0)

    print("Save locations...")
    # Save a file for each location generating pivot table
    # Exclude the last register because is noisy
    to_date = datetime.now() - timedelta(days=14)
    df_fit_list = []
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

        if location["country"] in df_cases_data["location"].values:

            df_fit = df_cases_data[df_cases_data["location"] == location["country"]].merge(df_location,
                                                                                           on=['location', 'date'],
                                                                                           how='outer')

            df_fit = df_fit.sort_values(['location', 'date', 'cases'])

            # Fit
            columns = [x for x in df_fit.columns.tolist() if x not in ["date", "location", "cases"]]
            for v in columns:
                df_fit[v] = df_fit["cases"] * (df_fit[v] / 100)
                df_fit[v] = df_fit[v].fillna(0.0)
                df_fit[v] = df_fit[v].apply(np.ceil).astype(int)

            df_fit['Unknown'] = df_fit.loc[:, columns].sum(axis=1)
            df_fit['Unknown'] = df_fit["cases"] - df_fit['Unknown']
            df_fit['Unknown'] = df_fit['Unknown'].fillna(0.0)
            df_fit['Unknown'] = df_fit['Unknown'].apply(np.ceil).astype(int)
            df_fit.loc[df_fit["Unknown"] < 0, "Unknown"] = 0

            df_fit.drop(columns=['cases'], inplace=True)

            if df_fit.size > 0:
                if df_fit.tail(1).iloc[0]['date'] > to_date:
                    # Remove the last register because probably is noise
                    df_fit = df_fit[:-1]

            df_fit.to_csv(
                f"../data/{location['country']}_fit.csv", index=False, quoting=csv.QUOTE_ALL, decimal=",")

            df_fit_list.append(df_fit)

    print("Generate Variants World data...")
    df_world = df.groupby(['date', 'variant']).agg({'num_sequences': 'sum'}).reset_index()

    dfb_world = df_world.groupby(['date']).agg(
        {'num_sequences': 'sum'}).rename(columns={"num_sequences": "num_sequences_total"}).reset_index()

    df_world = pd.merge(df_world, dfb_world, on=['date'])

    df_world["perc_sequences"] = (df_world["num_sequences"] / df_world["num_sequences_total"]) * 100

    df_world = df_world.drop(df_world[df_world.perc_sequences.isnull()].index)

    df_world = df_world.sort_values(['date', 'variant'])

    df_world_pivoted = df_world.pivot(index=["date"], columns=["variant"],
                                      values="perc_sequences").reset_index()
    df_world_pivoted["location"] = "World"
    df_world_pivoted.insert(0, "location", df_world_pivoted.pop("location"))
    df_world_pivoted = df_world_pivoted.fillna(0)

    print("Save World.csv...")
    df_world_pivoted[:-1].to_csv("../data/World.csv", index=False, quoting=csv.QUOTE_ALL, decimal=",")

    # Generate World Fit
    print("Create location list...")
    df_world_fit = pd.concat(df_fit_list)

    # Get columns to melt
    df_world_fit = df_world_fit.melt(id_vars=["location", "date"], var_name="variant", value_name="cases")

    # Clear zeroes
    df_world_fit = df_world_fit[df_world_fit.cases != 0]
    df_world_fit = df_world_fit.drop(df_world_fit[df_world_fit.cases.isnull()].index)

    df_world_fit = df_world_fit.groupby(['date', 'variant']).agg({'cases': 'sum'}).reset_index()

    df_world_fit = df_world_fit.sort_values(['date', 'variant'])

    df_world_fit_pivoted = df_world_fit.pivot(index=["date"], columns=["variant"], values="cases").reset_index()
    df_world_fit_pivoted["location"] = "World"
    df_world_fit_pivoted.insert(0, "location", df_world_fit_pivoted.pop("location"))
    df_world_fit_pivoted = df_world_fit_pivoted.fillna(0)

    df_world_fit_pivoted[:-1].to_csv("../data/World_fit.csv", index=False, quoting=csv.QUOTE_ALL, decimal=",")

    print("Generate update file...")
    update_dict = {
        "last_update_utc": datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S.%f%z')
    }
    with open('../data/update.json', 'w') as upd_file:
        json.dump(update_dict, upd_file, indent=4)

    print("End.")


main()
