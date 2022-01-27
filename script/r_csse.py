import pandas as pd


def main():
    df_csse_locations = pd.read_csv(
        "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/" +
        "csse_covid_19_data/UID_ISO_FIPS_LookUp_Table.csv")

    iso_list = set(df_csse_locations["iso3"].tolist())

    # Get the sub-regions that are actually countries
    subregion_to_region_iso = []
    for location in iso_list:
        df_location = df_csse_locations[
            (df_csse_locations["iso3"] == location) &
            (~df_csse_locations["Province_State"].isnull()) &
            (df_csse_locations["Province_State"] != df_csse_locations["Country_Region"])
            ]
        if len(df_location) == 1:
            subregion_to_region_iso.append(location)

    df_csse_loc_to_separate = df_csse_locations[df_csse_locations["iso3"].isin(subregion_to_region_iso)]
    subregion_to_region = sorted(df_csse_loc_to_separate["Province_State"].tolist())

    print("\n* CSSE Subregion to Region list")
    print("  " + "\n  ".join(subregion_to_region))

    # Countries that are not countries

    no_countries = df_csse_locations[df_csse_locations["iso3"].isnull()]["Country_Region"].tolist()

    print("\n* CSSE Countries that are not countries")
    print("  " + "\n  ".join(no_countries))


if __name__ == '__main__':
    main()
