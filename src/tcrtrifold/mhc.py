import requests
import polars as pl
import threading
from difflib import SequenceMatcher
from Bio import SeqIO
from pathlib import Path
import itertools
import warnings
import tempfile
import os

_tmp_dir = tempfile.TemporaryDirectory(prefix="tcrtrifold_")
HLA_UNIPROT_CACHE = Path(_tmp_dir.name) / "hla_uniprot_cache.csv"
tmp_df = pl.DataFrame({name: [] for name in ["name", "sequence", "accession"]})
tmp_df.write_csv(Path(_tmp_dir.name) / "hla_uniprot_cache.csv")


class MHCCodeConverter:
    def __init__(self):
        pass

    def get_sequence(self, official_name):
        pass


class MHCSequenceConverter:
    def __init__(self):
        pass

    def get_mhc_allele(self, sequence, thresh=100, chain=None):
        pass


# https://www.uniprot.org/uniprotkb/P01911/entry
# https://www.uniprot.org/uniprotkb/P79483/entry
# https://www.uniprot.org/uniprotkb/Q30154/entry
DRB_EXTRACELLULAR_TOPOLOGICAL_DOMAIN = slice(29, 227)

# https://www.uniprot.org/uniprotkb/P01903/entry#sequences
DRA_SEQ = "MAISGVPVLGFFIIAVLMSAQESWAIKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYTPITNVPPEVTVLTNSPVELREPNVLICFIDKFTPPVVNVTWLRNGKPVTTGVSETVFLPREDHLFRKFHYLPFLPSTEDVYDCRVEHWGLDEPLLKHWEFDAPSPLPETTENVVCALGLTVGLVGIIIGTIFIIKGLRKSNAAERRGPL"
DRA_NAME = "DRA*01:02"

# https://www.uniprot.org/uniprotkb/P01903/entry
DRA_EXTRACELLULAR_TOPOLOGICAL_SEQ = DRA_SEQ[26:216]

# https://www.uniprot.org/uniprotkb/P01906/entry
DQA2_EXTRACELLULAR_TOPOLOGICAL_DOMAIN = slice(23, 217)

# https://www.uniprot.org/uniprotkb/P01909/entry
DQA1_EXTRACELLULAR_TOPOLOGICAL_DOMAIN = slice(23, 216)

# https://www.uniprot.org/uniprotkb/P01920/entry
DQB1_EXTRACELLULAR_TOPOLOGICAL_DOMAIN = slice(32, 230)

# https://www.uniprot.org/uniprotkb/P05538/entry
DQB2_EXTRACELLULAR_TOPOLOGICAL_DOMAIN = slice(33, 229)

# https://www.uniprot.org/uniprotkb/P04439/entry
A_EXTRACELLULAR_TOPOLOGICAL_DOMAIN = slice(24, 308)

# https://www.uniprot.org/uniprotkb/P01889/entry
B_EXTRACELLULAR_TOPOLOGICAL_DOMAIN = slice(24, 309)

# https://www.uniprot.org/uniprotkb/P10321/entry
C_EXTRACELLULAR_TOPOLOGICAL_DOMAIN = slice(24, 308)

# https://www.uniprot.org/uniprotkb/P13747/entry
E_EXTRACELLULAR_TOPOLOGICAL_DOMAIN = slice(21, 305)

# https://www.uniprot.org/uniprotkb/P04440/entry#subcellular_location
DPB1_EXTRACELLULAR_TOPOLOGICAL_DOMAIN = slice(29, 225)

# https://www.uniprot.org/uniprotkb/P20036/entry
DPA1_EXTRACELLULAR_TOPOLOGICAL_DOMAIN = slice(28, 222)

DQA_FOR = {
    "DQB1*06:02": "DQA1*01:02",
    "DQB1*02:01": "DQA1*05:01",
    "DQB1*02:02": "DQA1*02:01",
    "DQB1*03:01": "DQA1*05:01",
    "DQB1*05:01": "DQA1*01:01",
    "DQB1*03:02": "DQA1*03:01",
    "DQB1*03:01": "DQA1*03:02",
    "DQB1*06:03": "DQA1*01:03",
    "DQB1*03:03": "DQA1*02:01",
    "DQB1*06:04": "DQA1*01:02",
    "DQB1*04:02": "DQA1*04:01",
    "DQB1*05:03": "DQA1*01:04",
    "DQB1*05:02": "DQA1*01:02",
}

# https://pmc.ncbi.nlm.nih.gov/articles/PMC2935290/
DPA_FOR = {
    "DPB1*04:01": "DPA1*01:03",
    "DPB1*01:01": "DPA1*02:01",
    "DPB1*02:01": "DPA1*01:03",
    "DPB1*04:02": "DPA1*03:01",
    "DPB1*05:01": "DPA1*02:01",
}

B2M_HUMAN_SEQ = "MSRSVALAVLALLSLSGLEAIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"


def _hla_slice_to_topological(official_name, sequence):
    if official_name.startswith("DRB"):
        return sequence[DRB_EXTRACELLULAR_TOPOLOGICAL_DOMAIN]
    elif official_name.startswith("DRA"):
        return DRA_EXTRACELLULAR_TOPOLOGICAL_SEQ
    elif official_name == "B2M":
        return sequence
    elif official_name.startswith("DQA1"):
        return sequence[DQA1_EXTRACELLULAR_TOPOLOGICAL_DOMAIN]
    elif official_name.startswith("DQA2"):
        return sequence[DQA2_EXTRACELLULAR_TOPOLOGICAL_DOMAIN]
    elif official_name.startswith("DQB1"):
        return sequence[DQB1_EXTRACELLULAR_TOPOLOGICAL_DOMAIN]
    elif official_name.startswith("DQB2"):
        return sequence[DQB2_EXTRACELLULAR_TOPOLOGICAL_DOMAIN]
    elif official_name.startswith("A"):
        return sequence[A_EXTRACELLULAR_TOPOLOGICAL_DOMAIN]
    elif official_name.startswith("B"):
        return sequence[B_EXTRACELLULAR_TOPOLOGICAL_DOMAIN]
    elif official_name.startswith("C"):
        return sequence[C_EXTRACELLULAR_TOPOLOGICAL_DOMAIN]
    elif official_name.startswith("E"):
        return sequence[E_EXTRACELLULAR_TOPOLOGICAL_DOMAIN]
    elif official_name.startswith("DPB1"):
        return sequence[DPB1_EXTRACELLULAR_TOPOLOGICAL_DOMAIN]
    elif official_name.startswith("DPA1"):
        return sequence[DPA1_EXTRACELLULAR_TOPOLOGICAL_DOMAIN]
    else:
        raise ValueError("No topological domain for sequence")


class HLACodeWebConverter(MHCCodeConverter):
    def __init__(self):

        self.cache = pl.read_csv(HLA_UNIPROT_CACHE)
        self.cache_lock = threading.Lock()

    def get_sequence(self, official_name, top_only=False):
        """
        https://hla.alleles.org/nomenclature/naming.html
        """
        if official_name == "B2M":
            return B2M_HUMAN_SEQ
        elif official_name.startswith("DRA"):
            sequence = DRA_SEQ
        else:
            with self.cache_lock:
                # search cache
                search = self.cache.filter(pl.col("name") == official_name).select(
                    "sequence"
                )
                if search.height == 1:
                    sequence = search.item()
                else:

                    params = {
                        "query": f"or(startsWith(name, '{official_name}:'), eq(name, '{official_name}'))",
                    }
                    # https://www.ebi.ac.uk/ipd/imgt/hla/alleles/
                    result = requests.get(
                        "https://www.ebi.ac.uk/cgi-bin/ipd/api/allele",
                        params=params,
                    )
                    result.raise_for_status()

                    json_result = result.json()

                    if json_result["meta"]["total"] == 0:
                        raise ValueError(
                            f"No allele found in IPD-IMGT/HLA for {official_name}"
                        )
                    # don't include unexpressed alleles
                    accessions = [
                        json_result["data"][i]["accession"]
                        for i in range(len(json_result["data"]))
                        if not json_result["data"][i]["name"].endswith("N")
                    ]
                    sequence = None

                    for accession in accessions:
                        result = requests.get(
                            f"https://www.ebi.ac.uk/cgi-bin/ipd/api/allele/{accession}"
                        )

                        result.raise_for_status()

                        json_result = result.json()
                        if json_result["status"] == "Public":
                            sequence = json_result["sequence"]["protein"]
                            break

                    if sequence is None:
                        raise ValueError(
                            f"No public sequence found for name {official_name}"
                        )

                    # cache the result
                    to_cache = pl.DataFrame(
                        {
                            "name": [official_name],
                            "sequence": [sequence],
                            "accession": [accession],
                        }
                    )
                    self.cache = pl.concat([self.cache, to_cache])
                    self.cache.write_csv(HLA_UNIPROT_CACHE)

        if top_only:
            sequence = _hla_slice_to_topological(official_name, sequence)

        return sequence


class HLASequenceDBConverter:
    """
    Converts an MHC sequence to its most likely HLA allele
    based on longest common substring
    """

    A_FNAMES = [
        "fasta/DQA1_prot.fasta",
        "fasta/DQA2_prot.fasta",
        "fasta/DRA_prot.fasta",
        "fasta/DPA1_prot.fasta",
    ]

    B_FNAMES = [
        "fasta/DRB1_prot.fasta",
        "fasta/DRB345_prot.fasta",
        "fasta/DQB1_prot.fasta",
        "fasta/DQB2_prot.fasta",
        "fasta/DPB1_prot.fasta",
    ]

    H_FNAMES = [
        "fasta/A_prot.fasta",
        "fasta/B_prot.fasta",
        "fasta/C_prot.fasta",
        "fasta/E_prot.fasta",
    ]

    STATUS = "Allele_status.txt"

    def __init__(self, db_repo_root):
        self.db_repo_root = Path(db_repo_root)

        self.alpha_seqs = self._load_fastas(self.A_FNAMES)
        self.beta_seqs = self._load_fastas(self.B_FNAMES)
        self.heavy_seqs = self._load_fastas(self.H_FNAMES)
        self.light_seqs = [("B2M", B2M_HUMAN_SEQ)]

        self._status_csv = pl.read_csv(
            self.db_repo_root / self.STATUS, comment_prefix="#"
        )

    def _get_status(self, name):
        search = self._status_csv.filter(pl.col("Allele") == name).select("Partial")
        if search.height != 1:
            raise ValueError(f"Expected 1 row in status csv for allele {name}")

        return search.item()

    def _load_fastas(self, fnames):
        ls = []
        for fname in fnames:
            with open(self.db_repo_root / fname) as f:
                fastas = list(SeqIO.parse(f, "fasta"))
                for fasta in fastas:
                    name, match_seq = fasta.description.split(" ")[1], str(fasta.seq)

                    # unexpressed allele
                    if name.endswith("N"):
                        continue

                    # decrease name to allele resolution
                    # name = ":".join(name.split(":")[:2])

                    ls.append((name, match_seq))

        return ls

    def get_mhc_allele(self, sequence, thresh=90, chain=None, top_only=False):
        """
        Returns the most likely HLA allele for a given sequence
        """
        # construct a list of candidate sequences
        potential_match = []
        if chain == "alpha":
            potential_match.extend(self.alpha_seqs)
        elif chain == "beta":
            potential_match.extend(self.beta_seqs)
        elif chain == "heavy":
            potential_match.extend(self.heavy_seqs)
        elif chain == "light":
            # there is only 1 chain possible
            return {
                "name": "B2M",
                "seq": B2M_HUMAN_SEQ,
                "max_resolution_name": "B2M",
                "sequence_status": "Full",
                "match_size": None,
            }
        else:
            potential_match.extend(self.alpha_seqs)
            potential_match.extend(self.beta_seqs)
            potential_match.extend(self.heavy_seqs)
            potential_match.extend(self.light_seqs)

        best_match = None
        best_seq = None
        best_name = None
        best_status = None

        for name, match_seq in potential_match:

            m = SequenceMatcher(
                None, sequence, match_seq, autojunk=False
            ).find_longest_match()

            if best_match is None or m.size > best_match.size:
                best_match = m
                best_name = name
                best_seq = match_seq
                best_status = self._get_status(name)

            elif m.size == best_match.size:
                # if we found an equally long match but with
                # "Full" status, prefer it
                new_status = self._get_status(name)
                if best_status == "Partial" and new_status == "Full":
                    best_match = m
                    best_name = name
                    best_seq = match_seq
                    best_status = new_status

        if best_match is None:
            raise ValueError(
                f"No matching allele found for sequence {sequence} with "
                f"threshold {thresh}"
            )

        # only return required resolution
        out_name = shorten_to_fullname(best_name)
        if top_only:
            out_seq = _hla_slice_to_topological(out_name, best_seq)

        else:
            out_seq = best_seq

        return {
            "name": out_name,
            "seq": out_seq,
            "max_resolution_name": best_name,
            "sequence_status": best_status,
            "match_size": best_match.size,
        }


def shorten_to_fullname(mhc_name):
    colon = mhc_name.count(":")
    if colon == 1:
        return mhc_name
    elif colon >= 2:
        return mhc_name.split(":")[0] + ":" + mhc_name.split(":")[1]
    else:
        return mhc_name


def is_fullname(mhc_name):
    return mhc_name.count(":") == 1


H2_I_HEAVY_DICT = {
    # https://www.ncbi.nlm.nih.gov/protein/1JUF_A?report=fasta
    # https://www.uniprot.org/uniprotkb/P01899/entry
    "H2-Db": "MGAMAPRTLLLLLAAALAPTQTRAGPHSMRYFETAVSRPGLEEPRYISVGYVDNKEFVRFDSDAENPRYEPRAPWMEQEGPEYWERETQKAKGQEQWFRVSLRNLLGYYNQSAGGSHTLQQMSGCDLGSDWRLLRGYLQFAYEGRDYIALNEDLKTWTAADMAAQITRRKWEQSGAAEHYKAYLEGECVEWLHRYLKNGNATLLRTDSPKAHVTHHPRSKGEVTLRCWALGFYPADITLTWQLNGEELTQDMELVETRPAGDGTFQKWASVVVPLGKEQNYTCRVYHEGLPEPLTLRWEPPPSTDSYMVIVAVLGVLGAMAIIGAVVAFVMKRRRNTGGKGGDYALAPGSQSSEMSLRDCKA",
    # https://www.ncbi.nlm.nih.gov/protein/3PAB_D?report=fasta
    # https://www.uniprot.org/uniprotkb/P01901/entry
    "H2-Kb": "MVPCTLLLLLAAALAPTQTRAGPHSLRYFVTAVSRPGLGEPRYMEVGYVDDTEFVRFDSDAENPRYEPRARWMEQEGPEYWERETQKAKGNEQSFRVDLRTLLGYYNQSKGGSHTIQVISGCEVGSDGRLLRGYQQYAYDGCDYIALNEDLKTWTAADMAALITKHKWEQAGEAERLRAYLEGTCVEWLRRYLKNGNATLLRTDSPKAHVTHHSRPEDKVTLRCWALGFYPADITLTWQLNGEELIQDMELVETRPAGDGTFQKWASVVVPLGKEQYYTCHVYHQGLPEPLTLRWEPPPSTVSNMATVAVLVVLGAAIVTGAVVAFVMKMRRRNTGGKGGDYALAPGSQTSDLSLPDCKVMVHDPHSLA",
    # https://www.ncbi.nlm.nih.gov/protein/5TRZ_C?report=fasta
    # https://www.uniprot.org/uniprotkb/P01902/entry
    "H2-Kd": "MAPCTLLLLLAAALAPTQTRAGPHSLRYFVTAVSRPGLGEPRFIAVGYVDDTQFVRFDSDADNPRFEPRAPWMEQEGPEYWEEQTQRAKSDEQWFRVSLRTAQRYYNQSKGGSHTFQRMFGCDVGSDWRLLRGYQQFAYDGRDYIALNEDLKTWTAADTAALITRRKWEQAGDAEYYRAYLEGECVEWLRRYLELGNETLLRTDSPKAHVTYHPRSQVDVTLRCWALGFYPADITLTWQLNGEDLTQDMELVETRPAGDGTFQKWAAVVVPLGKEQNYTCHVHHKGLPEPLTLRWKLPPSTVSNTVIIAVLVVLGAAIVTGAVVAFVMKMRRNTGGKGVNYALAPGSQTSDLSLPDGKVMVHDPHSLA",
    # https://www.ncbi.nlm.nih.gov/protein/CAA24128.1?report=fasta
    # https://www.uniprot.org/uniprotkb/P01897/entry
    "H2-Ld": "MGAMAPRTLLLLLAAALAPTQTRAGPHSMRYFETAVSRPGLGEPRYISVGYVDNKEFVRFDSDAENPRYEPQAPWMEQEGPEYWERITQIAKGQEQWFRVNLRTLLGYYNQSAGGTHTLQWMYGCDVGSDGRLLRGYEQFAYDGCDYIALNEDLKTWTAADMAAQITRRKWEQAGAAEYYRAYLEGECVEWLHRYLKNGNATLLRTDSPKAHVTHHPRSKGEVTLRCWALGFYPADITLTWQLNGEELTQDMELVETRPAGDGTFQKWASVVVPLGKEQNYTCRVYHEGLPEPLTLRWEPPPSTDSYMVIVAVLGVLGAMAIIGAVVAFVMKRRRNTGGKGGDYALAPGSQSSEMSLRDCKA",
    # https://www.uniprot.org/uniprotkb/P01900/entry
    "H2-Dd": "MGAMAPRTLLLLLAAALGPTQTRAGSHSLRYFVTAVSRPGFGEPRYMEVGYVDNTEFVRFDSDAENPRYEPRARWIEQEGPEYWERETRRAKGNEQSFRVDLRTALRYYNQSAGGSHTLQWMAGCDVESDGRLLRGYWQFAYDGCDYIALNEDLKTWTAADMAAQITRRKWEQAGAAERDRAYLEGECVEWLRRYLKNGNATLLRTDPPKAHVTHHRRPEGDVTLRCWALGFYPADITLTWQLNGEELTQEMELVETRPAGDGTFQKWASVVVPLGKEQKYTCHVEHEGLPEPLTLRWGKEEPPSSTKTNTVIIAVPVVLGAVVILGAVMAFVMKRRRNTGGKGGDYALAPGSQSSDMSL",
}


H2_I_HEAVY_TOP_DICT = {
    # https://www.uniprot.org/uniprotkb/P01899/entry
    "H2-Db": slice(24, 309),
    # https://www.uniprot.org/uniprotkb/P01901/entry
    "H2-Kb": slice(21, 305),
    # https://www.uniprot.org/uniprotkb/P01902/entry
    "H2-Kd": slice(21, 305),
    # https://www.uniprot.org/uniprotkb/P01897/entry
    "H2-Ld": slice(24, 309),
    # https://www.uniprot.org/uniprotkb/P01900/entry
    "H2-Dd": slice(24, 311),
}

H2_I_LIGHT_DICT = {
    # https://www.ncbi.nlm.nih.gov/protein/3P9M_B?report=fasta
    "B2M": "IQKTPQIQVYSRHPPENGKPNILNCYVTQFHPPHIEIQMLKNGKKIPKVEMSDMSFSKDWSFYILAHTEFTPTETDTYACRVKHDSMAEPKTVYWDRDM"
}

H2_I_LIGHT_TOP_DICT = {
    # https://www.uniprot.org/uniprotkb/P01887/entry
    "B2M": slice(None, None, None)
}
"MPRSRALILGVLALTTMLSLCGGEDDIEADHVGTYGISVYQSPGDIGQYTFEFDGDELFYVDLDKKETVWMLPEFGQLASFDPQGGLQNIAVVKHNLGVLTKRSNSTPATNEAPQATVFPKSPVLLGQPNTLICFVDNIFPPVINITWLRNSKSVADGVYETSFFVNRDYSFHKLSYLTFIPSDDDIYDCKVEHWGLEEPVLKHWEPEIPAPMSELTETVVCALGLSVGLVGIVVGTIFIIQGLRSGGTSRHPGPL"
H2_II_ALPHA_DICT = {
    # https://www.ncbi.nlm.nih.gov/protein/3MBE_E?report=fasta
    # https://www.uniprot.org/uniprotkb/P04228/entry
    "H2-IAg7": "MPCSRALILGVLALNTMLSLCGGEDDIEADHVGFYGTTVYQSPGDIGQYTHEFDGDELFYVDLDKKKTVWRLPEFGQLILFEPQGGLQNIAAEKHNLGILTKRSNFTPATNEAPQATVFPKSPVLLGQPNTLICFVDNIFPPVINITWLRNSKSVTDGVYETSFLVNRDHSFHKLSYLTFIPSDDDIYDCKVEHWGLEEPVLKHWEPEIPAPMSELTETVVCALGLSVGLVGIVVGTIFIIQGLRSGGTSRHPGPL",
    # https://www.ncbi.nlm.nih.gov/protein/AAA39548.1?report=fasta
    # https://www.uniprot.org/uniprotkb/P14434/entry
    "H2-IAb": "MPRSRALILGVLALTTMLSLCGGEDDIEADHVGTYGISVYQSPGDIGQYTFEFDGDELFYVDLDKKETVWMLPEFGQLASFDPQGGLQNIAVVKHNLGVLTKRSNSTPATNEAPQATVFPKSPVLLGQPNTLICFVDNIFPPVINITWLRNSKSVADGVYETSFFVNRDYSFHKLSYLTFIPSDDDIYDCKVEHWGLEEPVLKHWEPEIPAPMSELTETVVCALGLSVGLVGIVVGTIFIIQGLRSGGTSRHPGPL",
    # https://www.rcsb.org/sequence/1IEA fasta file
    # matches https://www.uniprot.org/uniprotkb/P01904/entry
    "H2-IEk": "MATIGALVLRFFFIAVLMSSQKSWAIKEEHTIIQAEFYLLPDKRGEFMFDFDGDEIFHVDIEKSETIWRLEEFAKFASFEAQGALANIAVDKANLDVMKERSNNTPDANVAPEVTVLSRSPVNLGEPNILICFIDKFSPPVVNVTWLRNGRPVTEGVSETVFLPRDDHLFRKFHYLTFLPSTDDFYDCEVDHWGLEEPLRKTWEFEEKTLLPETKENVMCALGLFVGLVGIVVGIILIMKGIKKRNVVERRQGAL",
}

H2_II_ALPHA_TOP_DICT = {
    # https://www.uniprot.org/uniprotkb/P04228/entry
    "H2-IAg7": slice(23, 218),
    # https://www.uniprot.org/uniprotkb/P14434/entry
    "H2-IAb": slice(23, 218),
    "H2-IEk": slice(25, 216),
}

H2_II_BETA_DICT = {
    # https://www.ncbi.nlm.nih.gov/protein/3MBE_F?report=fasta
    # Poor match in Uniprot
    "H2-IAg7": "GSGSGSGDSERHFVHQFKGECYFTNGTQRIRLVTRYIYNREEYLRFDSDVGEYRAVTELGRHSAEYYNKQYLERTRAELDTACRHNYEETEVPTSLRRLEQPNVAISLSRTEALNHHNTLVCSVTDFYPAKIKVRWFRNGQEETVGVSSTQLIRNGDWTFQVLVMLEMTPHQGEVYTCHVEHPSLKSPITVEWSSADLVPR",
    # https://www.ncbi.nlm.nih.gov/protein/AAA39548.1?report=fasta
    # https://www.uniprot.org/uniprotkb/P01921/entry
    "H2-IAb": "MALQIPSLLLSAAVVVLMVLSSPRTEGGNSERHFVVQFKGECYYTNGTQRIRLVTRYIYNREEYVRYDSDVGEYRAVTELGRPDAEYWNSQPEILERTRAEVDTACRHNYEGPETSTSLRRLEQPNVAISLSRTEALNHHNTLVCSVTDFYPAKIKVRWFRNGQEETVGVSSTQLIRNGDWTFQVLVMLEMTPHQGEVYTCHVEHPSLKSPITVEWRAQSESARSKMLSGIGGCVLGVIFLGLGLFIRHRSQKGPRGPPPAGLLQ",
    # https://www.uniprot.org/uniprotkb/P04230/entry
    # alias of Eb1
    "H2-IEk": "MVWLPRVPCVAAVILLLTVLSPPMALVRDSRPWFLEYCKSECHFYNGTQRVRLLERYFYNLEENLRFDSDVGEFHAVTELGRPDAENWNSQPEFLEQKRAEVDTVCRHNYEISDKFLVRRRVEPTVTVYPTKTQPLEHHNLLVCSVSDFYPGNIEVRWFRNGKEEKTGIVSTGLVRNGDWTFQTLVMLETVPQSGEVYTCQVEHPSLTDPVTVEWKAQSTSAQNKMLSGVGGFVLGLLFLGAGLFIYFRNQKGQSGLQPTGLLS",
}


H2_II_BETA_TOP_DICT = {
    # This indicates we don't know the topological domain/ don't have a full sequence
    "H2-IAg7": None,
    "H2-IAb": slice(27, 226),
    "H2-IEk": slice(26, 225),
}


def _h2_slice_to_topological(official_name, sequence, chain):
    if official_name == "B2M":
        return sequence
    elif official_name.startswith("H2-I"):
        if chain == "alpha":
            return sequence[H2_II_ALPHA_TOP_DICT[official_name]]
        elif chain == "beta":
            return sequence[H2_II_BETA_TOP_DICT[official_name]]
    else:
        return sequence[H2_I_HEAVY_TOP_DICT[official_name]]


class H2SequenceDictConverter:

    def __init__(self):
        pass

    def get_mhc_allele(self, sequence, thresh=90, chain=None, top_only=False):

        if chain is None:
            raise ValueError("Chain must be specified")
        # construct a list of candidate sequences
        potential_match = []
        if chain == "alpha":
            potential_match.extend(H2_II_ALPHA_DICT.items())
        elif chain == "beta":
            potential_match.extend(H2_II_BETA_DICT.items())
        elif chain == "heavy":
            potential_match.extend(H2_I_HEAVY_DICT.items())
        elif chain == "light":
            return {
                "name": "B2M",
                "seq": H2_I_LIGHT_DICT["B2M"],
                "max_resolution_name": "B2M",
                "sequence_status": "Full",
                "match_size": None,
            }
        else:
            raise ValueError(f"Unrecognized chain {chain}")

        best_match = None
        best_seq = None
        best_name = None

        for name, match_seq in potential_match:

            m = SequenceMatcher(
                None, sequence, match_seq, autojunk=False
            ).find_longest_match()

            if best_match is None or m.size > best_match.size:
                best_match = m
                best_name = name
                best_seq = match_seq

        if best_match is None:
            raise ValueError(
                f"No matching allele found for sequence {sequence} with "
                f"threshold {thresh}"
            )

        if top_only:
            out_seq = _h2_slice_to_topological(best_name, best_seq, chain)

        else:
            out_seq = best_seq

        return {
            "name": best_name,
            "seq": out_seq,
            "max_resolution_name": best_name,
            "sequence_status": "Full",
            "match_size": best_match.size,
        }


class H2CodeDictConverter(MHCCodeConverter):
    # https://biology.stackexchange.com/questions/76589/mouse-h2-allele-sequences
    # https://www.informatics.jax.org/
    def get_sequence(self, official_name, chain, top_only=False):
        if chain == "light":
            d = H2_I_LIGHT_DICT
            top = H2_I_LIGHT_TOP_DICT
        elif chain == "heavy":
            d = H2_I_HEAVY_DICT
            top = H2_I_HEAVY_TOP_DICT
        elif chain == "alpha":
            d = H2_II_ALPHA_DICT
            top = H2_II_ALPHA_TOP_DICT
        elif chain == "beta":
            d = H2_II_BETA_DICT
            top = H2_II_BETA_TOP_DICT
        else:
            raise ValueError(f"Unrecognized chain {chain}")

        if official_name not in d.keys():
            raise ValueError(f"Name {official_name} not found for chain {chain}")

        if top_only:
            if top[official_name] is None:
                warnings.warn(f"No topological domain known for {official_name}")
                return d[official_name]
            return d[official_name][top[official_name]]

        return d[official_name]
