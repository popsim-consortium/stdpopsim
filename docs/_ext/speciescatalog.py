import csv
import pathlib

from docutils import nodes
from docutils.parsers.rst import Directive

import stdpopsim


def species_summary_table(species):
    table = nodes.table()
    tgroup = nodes.tgroup(cols=2)
    colspec = nodes.colspec(colwidth=1)
    tgroup.append(colspec)
    colspec = nodes.colspec(colwidth=1)
    tgroup.append(colspec)
    colspec = nodes.colspec(colwidth=1)
    tgroup.append(colspec)

    table += tgroup

    data = [
        ("id", species.id, ""),
        ("name", species.name, ""),
        ("generation_time", species.generation_time, "TODO: Notes for generation_time"),
        ("population_size", species.population_size, "TODO: Notes for population_size"),
    ]

    rows = []
    for row_data in data:
        row = nodes.row()
        rows.append(row)
        for entry_data in row_data:
            entry = nodes.entry()
            entry += nodes.paragraph(text=entry_data)
            row += entry

    tbody = nodes.tbody()
    tbody.extend(rows)
    tgroup += tbody

    return table


def model_table(model):
    table = nodes.table()
    tgroup = nodes.tgroup(cols=2)

    colspec = nodes.colspec(colwidth=1)
    tgroup.append(colspec)
    colspec = nodes.colspec(colwidth=1)
    tgroup.append(colspec)

    table += tgroup

    rows = []

    row = nodes.row()
    rows.append(row)
    entry = nodes.entry()
    entry += nodes.paragraph(text="kind")
    row += entry
    entry = nodes.entry()
    entry += nodes.paragraph(text=model.kind)
    row += entry

    row = nodes.row()
    rows.append(row)
    entry = nodes.entry()
    entry += nodes.paragraph(text="name")
    row += entry
    entry = nodes.entry()
    entry += nodes.paragraph(text=model.name)
    row += entry

    row = nodes.row()
    rows.append(row)
    entry = nodes.entry()
    entry += nodes.paragraph(text="num_populations")
    row += entry
    entry = nodes.entry()
    entry += nodes.paragraph(text=model.num_populations)
    row += entry

    row = nodes.row()
    rows.append(row)
    entry = nodes.entry()
    entry += nodes.paragraph(text="doi")
    row += entry
    entry = nodes.entry()
    para = nodes.paragraph()
    para += nodes.reference(internal=False, refuri=model.doi, text=model.doi)
    entry += para
    row += entry

    tbody = nodes.tbody()
    tbody.extend(rows)
    tgroup += tbody

    return table


def citation_list(citable):
    bullet_list = nodes.bullet_list()
    for citation in citable.citations:
        list_item = nodes.list_item()
        list_item += nodes.paragraph(text=citation)
        bullet_list += list_item
    return bullet_list


def model_parameter_table(model):
    path = pathlib.Path(f"parameter_tables/{model.species.name}/{model.name}.csv")
    if not path.exists():
        return None
    with open(path) as csv_file:
        reader = csv.reader(csv_file)
        data = list(reader)

    table = nodes.table()
    tgroup = nodes.tgroup(cols=3)
    colspec = nodes.colspec(colwidth=1)
    tgroup.append(colspec)
    colspec = nodes.colspec(colwidth=1)
    tgroup.append(colspec)
    colspec = nodes.colspec(colwidth=1)
    tgroup.append(colspec)

    table += tgroup

    thead = nodes.thead()
    tgroup += thead
    row = nodes.row()
    entry = nodes.entry()
    entry += nodes.paragraph(text="Parameter Type (units)")
    row += entry
    entry = nodes.entry()
    entry += nodes.paragraph(text="Value")
    row += entry
    entry = nodes.entry()
    entry += nodes.paragraph(text="Description")
    row += entry
    thead.append(row)

    rows = []
    for row_data in data:
        row = nodes.row()
        rows.append(row)
        for entry_data in row_data:
            entry = nodes.entry()
            entry += nodes.paragraph(text=entry_data)
            row += entry

    tbody = nodes.tbody()
    tbody.extend(rows)
    tgroup += tbody

    return table


def population_list(model):
    pop_list = nodes.enumerated_list(start=0)
    for population in model.populations:
        list_item = nodes.list_item()
        # TODO add the population name in Strong here somehow.
        list_item += nodes.paragraph(text=population.description)
        pop_list += list_item
    return pop_list


def genetic_map_section(species, genetic_map):
    map_id = f"sec_catalog_{species.name}_genetic_maps_{genetic_map.name}"
    section = nodes.section(ids=[map_id])
    section += nodes.title(text=genetic_map.name)
    section += nodes.paragraph(text=genetic_map.description)
    return section

def genetic_maps_table(species):
    table = nodes.table()
    tgroup = nodes.tgroup(cols=2)
    colspec = nodes.colspec(colwidth=1)

    tgroup.append(colspec)
    colspec = nodes.colspec(colwidth=1)
    tgroup.append(colspec)
    colspec = nodes.colspec(colwidth=1)
    tgroup.append(colspec)

    table += tgroup

    thead = nodes.thead()
    tgroup += thead
    row = nodes.row()
    entry = nodes.entry()
    entry += nodes.paragraph(text="Name")
    row += entry
    entry = nodes.entry()
    entry += nodes.paragraph(text="Year")
    row += entry
    entry = nodes.entry()
    entry += nodes.paragraph(text="Reference")
    row += entry

    thead.append(row)

    rows = []
    for genetic_map in species.genetic_maps:
        row = nodes.row()
        rows.append(row)

        # map_id = f"sec_catalog_{species.name}_genetic_maps_{genetic_map.name}"
        # TODO figure out how to make a link here to the corresponding detail
        # section.
        entry = nodes.entry()
        entry += nodes.paragraph(text=genetic_map.name)
        row += entry

        entry = nodes.entry()
        entry += nodes.paragraph(text=genetic_map.year)
        row += entry

        entry = nodes.entry()
        para = nodes.paragraph()
        para += nodes.reference(
            internal=False, refuri=genetic_map.doi, text=genetic_map.doi)
        entry += para
        row += entry

    tbody = nodes.tbody()
    tbody.extend(rows)
    tgroup += tbody

    return table


def chromosomes_table(species):

    table = nodes.table()
    tgroup = nodes.tgroup(cols=2)
    colspec = nodes.colspec(colwidth=1)
    tgroup.append(colspec)
    colspec = nodes.colspec(colwidth=1)
    tgroup.append(colspec)

    table += tgroup

    thead = nodes.thead()
    tgroup += thead
    row = nodes.row()
    entry = nodes.entry()
    entry += nodes.paragraph(text="Name")
    row += entry

    entry = nodes.entry()
    entry += nodes.paragraph(text="Length")
    row += entry

    thead.append(row)

    rows = []
    for chrom in species.genome.chromosomes:
        row = nodes.row()
        entry = nodes.entry()
        entry += nodes.paragraph(text=chrom.name)
        row += entry

        entry = nodes.entry()
        entry += nodes.paragraph(text="{:.4G}".format(chrom.length))
        row += entry

        # TODO add mutation/recombination rate.
        rows.append(row)
    tbody = nodes.tbody()
    tbody.extend(rows)
    tgroup += tbody

    return table


def model_section(species, model):
    section = nodes.section(ids=[f"sec_catalog_{species.name}_models_{model.name}"])
    section += nodes.title(text=model.name)
    section += nodes.paragraph(text=model.description)
    section += nodes.rubric(text="Details")
    section += model_table(model)
    section += nodes.rubric(text="Populations")
    section += population_list(model)
    section += nodes.rubric(text="Citations")
    section += citation_list(model)
    section += nodes.rubric(text="Model parameters")
    section += model_parameter_table(model)
    return section


class SpeciesCatalog(Directive):

    required_arguments = 1

    def run(self):
        species = stdpopsim.get_species(self.arguments[0])
        sid = f"sec_catalog_{species.name}"
        section = nodes.section(ids=[sid])
        # TODO figure out out to make cross references to these IDs work within sphinx.
        # Probably need to create a 'label' node?
        # section += nodes.label(ids=[sid])
        section += nodes.title(text=species.name)
        section += species_summary_table(species)

        genome_section = nodes.section(ids=[f"sec_catalog_{species.name}_genome"])
        genome_section += nodes.title(text="Genome")
        genome_section += chromosomes_table(species)
        section += genome_section

        maps_section = nodes.section(ids=[f"sec_catalog_{species.name}_genetic_maps"])
        maps_section += nodes.title(text="Genetic Maps")
        maps_section += genetic_maps_table(species)
        for gmap in species.genetic_maps:
            maps_section += genetic_map_section(species, gmap)
        section += maps_section

        models_section = nodes.section(ids=[f"sec_catalog_{species.name}_models"])
        models_section += nodes.title(text="Models")
        # TODO add a table summarising the models with links to the detailed
        # descriptions
        for model in species.models:
            models_section += model_section(species, model)
        section += models_section

        return [section]


def setup(app):
    app.add_directive("speciescatalog", SpeciesCatalog)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
