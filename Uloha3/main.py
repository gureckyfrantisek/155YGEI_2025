from lines_to_graph2 import *
from pathlib import Path

# Implement Dijskra's algo here

def main():
    # First translate all line data to graphs
    data_dir = Path(__file__).parent / "data"
    graphs = []

    for file in data_dir.iterdir():
        graphs.append(convertToGraph(file))

if __name__ == "__main__":
    main()
