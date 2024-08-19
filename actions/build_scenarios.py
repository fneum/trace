#!/usr/bin/env python
# coding: utf-8

# In[14]:

from itertools import product
import pandas as pd
import yaml

# In[8]:

with open("config/config.default.yaml") as f:
    config = yaml.safe_load(f)
exporters = list(config["regions"].keys())


# In[9]:


importers = {
    "shipping": [
        "T0557",  # South Hook LNG terminal, UK
        "T0492",  # Gate LNG terminal, Netherlands
        "T0498",  # Swinoujscie LNG terminal, Poland
        "T0462",  # Revithoussa LNG terminal, Greece
        "T0466",  # Adriatic LNG terminal, Italy
        "T0522",  # Barcelona LNG terminal, Spain
        "T0500",  # Sines LNG terminal, Portugal
    ],
    "pipeline-cable": [
        "EUSW",
        "EUS",
        "EUSE",
        "EUE",
    ]
}


# In[10]:


escs = [
    "hvdc-to-elec",
    "pipeline-h2",
    "shipping-ftfuel",
    "shipping-lch4",
    "shipping-lh2",
    "shipping-lnh3",
    "shipping-meoh",
    "shipping-steel",
    "shipping-hbi",
]


# In[11]:


forbidden = {
    "shipping": [
        "KZ",
        "UA",
        "AZ",
        "TM",
        "UZ",
        "MN",
        "AF",
    ],
    "pipeline-cable": [
        "AR-South",
        "AR-North",
        "AU-West",
        "AU-East",
        "BR-Northeast",
        "BR-Southeast",
        "CA-East",
        "CL",
        "NA",
        "US-Northeast",
        "US-Southeast",
        "US-Northwest",
        "US-South",
        "US-Southwest",
        "US-Alaska",
        "MX",
        "ER",
        "PK",
        "VN",
        "MM",
        "TH",
        "LK",
        "ZA",
        "KE",
        "TZ",
        "NG",
        "AO",
        "ET",
        "MZ",
        "SN",
        "UY",
        "CO",
        "PE",
        "IN-Northwest",
        "IN-South",
        "IN-East",
        "SO",
        "MG",
        "DJ",
        "BO",
        "VE",
    ],
}


# In[29]:


forbidden_combinations = {
    "EUNW": set(exporters) - {"IS"},
    "EUSW": [
        "CN-Northeast",
        "CN-Southeast",
        "CN-West",
        "EG",
        "KZ",
        "LY",
        "OM",
        "SA",
        "TR",
        "UA",
        "GE",
        "AZ",
        "AE",
        "QA",
        "KW",
        "SY",
        "TM",
        "UZ",
        "MN",
        "IS",
        "IL",
        "SD",
        "YE",
        "JO",
        "IR",
        "IQ",
        "AF",
    ],
    "EUS": [
        "CN-Northeast",
        "CN-Southeast",
        "CN-West",
        "KZ",
        "EG",
        "MA",
        "EH",
        "OM",
        "SA",
        "TR",
        "UA",
        "GE",
        "AZ",
        "QA",
        "AE",
        "KW",
        "TM",
        "UZ",
        "MN",
        "MR",
        "IS",
        "IL",
        "SD",
        "SY",
        "YE",
        "JO",
        "IR",
        "IQ",
        "AF",
    ],
    "EUSE": [
        "DZ",
        "MA",
        "EH",
        "TN",
        "MR",
        "IS",
        "UA",
    ],
    "EUE": [
        "DZ",
        "EG",
        "LY",
        "MA",
        "EH",
        "OM",
        "SA",
        "TN",
        "TR",
        "QA",
        "AE",
        "KW",
        "MR",
        "IL",
        "SD",
        "SY",
        "YE",
        "JO",
        "IQ",
        "IS",
    ],
}


# In[30]:


# years = [2030, 2050]
years = [2030]
# scenarios = ["default", "low-electrolysis-50", "cavern-h2-storage"]
scenarios = ["default"]
scenario_list = []
for year, scenario, esc, exporter in product(years, scenarios, escs, exporters):

    if "shipping" in esc:
        forbidden_esc_exporters = forbidden["shipping"]
        esc_importers = importers["shipping"]
    else:
        forbidden_esc_exporters = forbidden["pipeline-cable"]
        esc_importers = importers["pipeline-cable"]
    if exporter in forbidden_esc_exporters:
        continue


    for esc_importer in esc_importers:

        if exporter in forbidden_combinations.get(esc_importer, []):
            continue
        
        scenario_list.append(pd.Series({
            "scenario": "default",
            "year": year,
            "esc": esc,
            "exporter": exporter,
            "importer": esc_importer,
        }))
scenario_list = pd.concat(scenario_list, axis=1).T
extras = pd.read_csv("scenarios/extra.csv")
scenario_list = pd.concat([scenario_list, extras], axis=0, ignore_index=True)
scenario_list.to_csv("scenarios/default.csv", index=False)

