import time
import requests
from bs4 import BeautifulSoup
import logging
import argparse
from pathlib import Path
from urllib.parse import urlencode


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def make_dir(location):
    Path(location).mkdir(parents=True, exist_ok=True)


def cleanup(folder):
    for file in folder.glob("*.fasta"):
        Path(file).unlink()
    Path(folder).rmdir()
    logging.info(
        f"Deleting folder: {folder.name}, because --cleanup was selected")


def create_library(folder: Path):
    library = folder.parent / "library"
    make_dir(library)
    with open(library / "library.fasta", 'w') as w:
        for file in folder.glob("*.fasta"):
            with open(file) as f:
                w.write(f.read())


def construct_url(segment, species, frame):
    base_url = "https://www.imgt.org/genedb/GENElect"
    params = {
        'query': f"{frame} {segment}",
        'species': species
    }
    return f"{base_url}?{urlencode(params)}"


def scrape(response):
    soup = BeautifulSoup(response.text, 'html.parser')
    paragraphs = soup.find_all('pre')
    for p in paragraphs:
        seq = p.text.strip()
        if ">" in seq:
            logging.info(f"Succeeded to retrieve sequences from IMGT")
            return seq


def write_sequence(name, folder, sequence):
    file = folder / f"{name}.fasta"
    logging.info(f"Writing sequences from {name} to {file.name}")
    with open(file, 'w') as f:
        f.write(str(sequence) + "\n")


def fetch_sequence(segment, folder, species, frame):
    url = construct_url(segment, species, frame)
    response = requests.get(url)
    if response.status_code == 200:
        sequence = scrape(response)
        if sequence:
            write_sequence(segment, folder, sequence)
        else:
            logging.info(f"No sequences found for {segment} of {species}.")
    else:
        logging.info(
            f"Failed to fetch data for {segment} of {species}. Status code: {response.status_code}")
    logging.info("Waiting 2 seconds to avoid overloading the IMGT server!")
    time.sleep(2)


def scrape_IMGT(species, imune_type, folder, frame):
    segments = {"TCR": ["TRBV", "TRBJ", "TRBD", "TRAV", "TRAJ", "TRDD", "TRDJ", "TRDV", "TRGV", "TRGJ"],
                "IG": ['IGHV', 'IGHD', 'IGHJ', 'IGKV', 'IGKJ', 'IGLV', 'IGLJ']}
    if not folder.exists():
        logging.info(f"Folder {folder.name} does not exist, creating it!")
        make_dir(folder)
    for segment in segments[imune_type]:
        logging.info(
            f"Retrieving sequences from IMGT for the {segment} of {species}")
        if not Path(folder / f"{segment}.fasta").exists():
            fetch_sequence(segment, folder, species, frame)
        else:
            logging.info(f"File {segment}.fasta already exists skipping!")
            time.sleep(2)


def to_capital(string: str):
    return string.capitalize()


def convert_frame(frame):
    options = {"all": "7.2", "in-frame": "7.5", "in-frame-gaps": "7.1"}
    return options[frame or "all"]


def argparser_setup():
    latin_names = [
        'Homo sapiens', 'Mus', 'Bos taurus',
        'Ovis aries', 'Danio rerio', 'Canis lupus familiaris',
        'Macaca mulatta', 'Tursiops truncatus', 'Oryctolagus cuniculus',
        'Heterocephalus glaber'
    ]
    parser = argparse.ArgumentParser(
        description='Scrape IMGT for TCR and IG segment sequences of a given species.')
    parser.add_argument('-S', "--species", type=to_capital, choices=latin_names, required=True,
                        help='Name of the species to scrape for (e.g., "Homo sapiens")')
    parser.add_argument('-T', '--type', type=str.upper, choices=[
                        'TCR', 'IG'], required=True, help='Choose between TCR (T-cell receptor) or IG (Immunoglobulin)')
    parser.add_argument('-O', '--output', type=str,
                        help='Output directory where the results will be saved.')
    parser.add_argument('-f', '--frame-selection', type=str, choices=['all', 'in-frame', 'in-frame-gaps'],
                        help='Select ORF frame analysis type: "all" for F+ORF+all P, "in-frame" for F+ORF+in-frame P , or "in-frame-gaps" for F+ORF+in-frame P with IMGT gaps.')
    parser.add_argument('--create-library', action='store_true',
                        help='Create a library from the IMGT files if specified.')
    parser.add_argument('--cleanup', action='store_true',
                        help='Clean up leftover IMGT files after processing.')

    args = parser.parse_args()
    return args


def main():
    cwd = Path.cwd()
    args = argparser_setup()
    logging.info(f"Selected species: {args.species} and type: {args.type}")
    folder_name = args.species.replace(" ", "_").lower()
    if args.output:
        folder_name = args.output
    folder = cwd / folder_name
    library = Path.cwd() / "library" / "library.fasta"
    if args.cleanup is not None and not Path.exists(library):
        scrape_IMGT(args.species, args.type, folder,
                    convert_frame(args.frame_selection))
        if not library.exists():
            if args.create_library:
                logging.info(f"Creating a library for generating files.")
                create_library(cwd / folder)
            if args.cleanup:
                cleanup(cwd / folder)
    else:
        logging.info(
            f"library.fasta already exists. If --cleanup was selected, please remove it if you want to downloads the files.")


if __name__ == '__main__':
    main()
