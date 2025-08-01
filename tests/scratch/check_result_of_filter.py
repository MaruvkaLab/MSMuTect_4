import os
from dataclasses import dataclass

import pandas as pd


@dataclass
class Case:
    MSI: str
    tumor_name: str
    normal_name: str

def relevant_line_numbers(fp: str):
    df = pd.read_csv(fp, delimiter="\t")
    noisy_col = df["Noisy Locus"]
    return noisy_col

def download_file(name: str, germline: str):
    os.system(f"gcloud storage cp gs://texas-bleeding/histogram_files/{name}.zip .")
    os.system(f"unzip -o {name}.zip -d {germline}")
    fname = f"{germline}/{name}.hist.tsv"
    return fname


def main():
    tumors_list = [
        "f04e96e5-e09e-4e79-ae5d-3f598ec41f29",
        "a19c9fef-711e-4885-ab68-fd57c57d278e",
        "bc368650-f168-487b-8eea-ef4817d3cb12",
        "f6c5cafd-6d29-430a-a3ca-d92c6eae456c",
        "178550ea-1cde-476e-ba2c-9f002193fd7f",
        "f24fbfd5-a9d9-45ab-9744-884a3599f47a",
        "6d82e4a5-1d2c-47b6-94a7-a43f22e0df83",
        "ee4f9de8-db99-4104-a120-7a75ec71d198",
        "2c3900fe-afc5-4be4-9949-ad758972898b",
        "7ba8ff74-9d29-4a4e-91b8-64cee9f223d9"
    ]

    normals_list = [

        "75f64f86-fab8-48be-97d6-bdcb11035aa6",
        "fa0aa386-66cd-408c-8ddd-9ddf0d9a5b59",
        "ad8c5132-383a-40d1-9691-51b898da0bb0",
        "d3d5f613-69a1-4d42-9447-f7126c168a8e",
        "7fabb56f-eb87-4a26-b5eb-24a9ca57fe19",
        "39d98988-6dc9-4dec-bf07-c7156ae46b59",
        "f05f2972-9a93-4aae-b118-ee05d7361756",
        "8082c449-619d-4870-8977-1d761c7a5b5d",
        "70a37d40-1b9b-438a-b6a7-a035fb23b1ee",
        "bcd43189-ba51-4b69-8d3e-a681108c399f"
    ]
    msi_list = ["MSI" for i in range(5)]+["MSS" for j in range(5)]
    cases = []
    for i in range(10):
        cases.append(Case(msi_list[i], tumors_list[i], normals_list[i]))

    for c in cases:
        normal_file = download_file(c.normal_name, "normal")
        tumor_file = download_file(c.tumor_name, "tumor")
        normal_noisy = relevant_line_numbers(normal_file)
        tumor_noisy = relevant_line_numbers(tumor_file)
        print(f"{c.MSI}: {((tumor_noisy+normal_noisy)>0).sum()}")
        os.remove(normal_file)
        os.remove(tumor_file)
        os.remove(f"{c.normal_name}.zip")
        os.remove(f"{c.tumor_name}.zip")


main()