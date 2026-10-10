"""Export a clean research-book scaffold without project history or data."""
from pathlib import Path
import argparse
import shutil


def export_project(destination, title):
    root = Path(__file__).resolve().parents[1]
    destination = Path(destination).expanduser().resolve()
    if destination == root or root in destination.parents:
        raise ValueError('Choose a destination outside the source repository.')
    if destination.exists() and any(destination.iterdir()):
        raise ValueError('Destination must be empty.')
    destination.mkdir(parents=True, exist_ok=True)
    names = ['_quarto.yml', '.gitignore', 'AGENTS.md', 'CONTRIBUTING.md',
             'requirements.txt', 'references.qmd', 'templates/investigation.qmd',
             '.github/workflows/book.yml', 'scripts/R/001-monte-carlo.R',
             'chapters/03-reproducible-example.qmd', 'scripts/new_project.py',
             'docs/new-project.md', 'scripts/python/README.md']
    for name in names:
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root / name, target)
    workflow = destination / '.github/workflows/book.yml'
    workflow.write_text(workflow.read_text().split('  publish:\n')[0].replace(
        'name: Build and publish research book', 'name: Build research book'))
    config = (destination / '_quarto.yml').read_text()
    import json
    config = '\n'.join(line for line in config.splitlines()
                       if not line.startswith('  site-url:')) + '\n'
    config = config.replace('"Alternative DAG Parameterizations"', json.dumps(title))
    (destination / '_quarto.yml').write_text(config)
    chapters = {
        '01-questions': 'Introduction and research questions {#sec-questions}',
        'related-work': 'Related work',
        '02-framework': 'Notation and theoretical framework {#sec-framework}',
        'methods': 'Methods and reproducibility',
        '04-findings': 'Research results {#sec-findings}',
        'discussion': 'Discussion',
        '05-open-questions': 'Conclusions and open questions',
        'research-log': 'Research log',
    }
    for name, heading in chapters.items():
        (destination / 'chapters' / (name + '.qmd')).write_text(
            '# ' + heading + '\n\nTo be developed. No research findings established yet.\n')
    (destination / 'index.qmd').write_text('# Overview {.unnumbered}\n\n' + title + '\n\n## Abstract\n\nResearch scope and results to be developed.\n')
    (destination / 'references.bib').write_text('% Add verified references here.\n')
    (destination / 'PROJECT_CONTEXT.md').write_text(
        '# Shared research context\n\nProject: ' + title +
        '\n\nRepository and website: not configured.\nScientific scope: to be defined.\nFindings: none.\n\nRead AGENTS.md and the research log at the start of each chat.\n')
    (destination / 'README.md').write_text(
        '# ' + title + '\n\nA Quarto research book with reusable R and Python investigations.\n\n'
        'Read CONTRIBUTING.md and PROJECT_CONTEXT.md before contributing.\n\n'
        '## Setup\n\nInstall Quarto and R 4.4.2.\n\n```sh\n'
        'Rscript -e \'install.packages(c("knitr", "rmarkdown", "renv"), repos="https://cloud.r-project.org")\'\n'
        'Rscript -e \'renv::init(bare=TRUE); renv::install(c("knitr", "rmarkdown")); renv::snapshot()\'\n'
        'quarto render\n```\n\n'
        'Commit the generated renv configuration and lockfile. Python dependencies belong in requirements.txt.\n\n'
        '## GitHub\n\nInitialize Git and connect a new repository. The exported workflow builds only. Add a project-specific publishing job and credentials if desired; no source-project deployment secrets or website targets are copied. Set book.site-url in _quarto.yml for hosting.\n\n'
        '## Layout\n\nchapters/: scientific narrative; scripts/: computations; materials/: source inventory; data/: input provenance; notes/: working reasoning; results/: generated outputs; references.bib: citations.\n\n'
        'Use one chat per investigation. Preserve substantive reasoning in chapters and the research log. Full chat history is not automatically shared.\n\n'
        'See docs/new-project.md to reuse this scaffold. No license has been selected.\n')
    for name in ['data/raw', 'materials/private', 'notes', 'results']:
        folder = destination / name
        folder.mkdir(parents=True, exist_ok=True)
        (folder / 'README.md').write_text('# ' + name + '\n\nRecord provenance and purpose before adding material.\n')
    print(destination)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination')
    parser.add_argument('--title', required=True)
    args = parser.parse_args()
    export_project(args.destination, args.title)
