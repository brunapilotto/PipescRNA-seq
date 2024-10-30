#!/usr/bin/env python3

import argparse
import json


def menu():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r", "--retinoblastoma_signature", help="Signature to update", required=True
    )
    parser.add_argument(
        "-i", "--imunno_signature", help="Immunologic sample signature", required=True
    )
    args = parser.parse_args()
    return args


def update_signature(retinoblastoma_signature: str, immunologic_signature: str) -> None:
    with open(retinoblastoma_signature) as retinoblastoma_file:
        base_file = json.load(retinoblastoma_file)
    with open(immunologic_signature) as immunologic_file:
        new_file = json.load(immunologic_file)

    for key, value in new_file.items():
        base_file[key] = {"geneSymbols": value}

    with open("updated_signature.json", "w") as updated_file:
        json.dump(base_file, updated_file, indent=4)


def main():
    args = menu()
    update_signature(args.retinoblastoma_signature, args.imunno_signature)


if __name__ == "__main__":
    main()
