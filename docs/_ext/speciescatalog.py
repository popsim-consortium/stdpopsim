"""
A sphinx extension that provides the "speciescatalog" directive. See the
Sphinx documentation for details on developing such extensions:

http://www.sphinx-doc.org/en/master/extdev/index.html#dev-extensions

This is based on the todo tutorial:
http://www.sphinx-doc.org/en/master/development/tutorials/todo.html
"""
import csv
import pathlib

from docutils import nodes
from sphinx.util.docutils import SphinxDirective

import stdpopsim


class SpeciesCatalogDirective(SphinxDirective):

    required_arguments = 1

    def get_target(self, tid):
        """
        Returns a target node for the specified ID.

        This took a lot of digging through the Docutils and Sphinx sources to
        figure out: it appears that including the 'names' argument to the
        target node and calling note_explicit_target are essential to allow
        us to cross reference these targets from the rest of the documentation.
        """
        target = nodes.target('', '', ids=[tid], names=[tid])
        self.state.document.note_explicit_target(target)
        return target

    def species_summary_table(self, species):
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
            ("generation_time", species.generation_time,
                "TODO: Notes for generation_time"),
            ("population_size", species.population_size,
                "TODO: Notes for population_size"),
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

    def model_table(self, model):
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
        entry += nodes.paragraph(text="id")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text=model.id)
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

        tbody = nodes.tbody()
        tbody.extend(rows)
        tgroup += tbody

        return table

    def citation_list(self, citable):
        bullet_list = nodes.bullet_list()
        for citation in citable.citations:
            list_item = nodes.list_item()

            para = nodes.paragraph(text=f"{citation.author}, {citation.year}. ")
            para += nodes.reference(
                internal=False, refuri=citation.doi, text=citation.doi)
            list_item += para
            bullet_list += list_item
        return bullet_list

    def model_parameter_table(self, species, model):
        path = pathlib.Path(f"parameter_tables/{species.id}/{model.id}.csv")
        if not path.exists():
            # There doesn't seem to be a better way to do this...
            print(
                f"FIXME: Skipping model parameters for {model.id} due to missing table")
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

    def population_list(self, model):
        pop_list = nodes.enumerated_list(start=0)
        for population in model.populations:
            list_item = nodes.list_item()
            # TODO add the population name in Strong here somehow.
            list_item += nodes.paragraph(text=population.description)
            pop_list += list_item
        return pop_list

    def genetic_map_section(self, species, genetic_map):
        # NOTE: Ids must be lowercase to work properly. When this was developed,
        # genetic maps didn't have standardised lowercase IDs. The .lower() here
        # should be removed in the future when this is sorted out.
        map_id = f"sec_catalog_{species.id}_genetic_maps_{genetic_map.name}".lower()
        target = self.get_target(map_id)
        section = nodes.section(ids=[map_id])
        section += nodes.title(text=genetic_map.name)
        section += nodes.paragraph(text=genetic_map.description)
        return [target, section]

    def genetic_maps_table(self, species):
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

            # map_id = f"sec_catalog_{species.id}_genetic_maps_{genetic_map.name}"
            # TODO figure out how to make a link here to the corresponding detail
            # section.
            entry = nodes.entry()
            entry += nodes.paragraph(text=genetic_map.name)
            row += entry

            entry = nodes.entry()
            entry += nodes.paragraph(text=genetic_map.citations[0].year)
            row += entry

            entry = nodes.entry()
            para = nodes.paragraph()
            doi = genetic_map.citations[0].doi
            para += nodes.reference(internal=False, refuri=doi, text=doi)
            entry += para
            row += entry

        tbody = nodes.tbody()
        tbody.extend(rows)
        tgroup += tbody

        return table

    def chromosomes_table(self, species):

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
            entry += nodes.paragraph(text=chrom.id)
            row += entry

            entry = nodes.entry()
            entry += nodes.paragraph(text="{:d}".format(chrom.length))
            row += entry

            # TODO add mutation/recombination rate.
            rows.append(row)
        tbody = nodes.tbody()
        tbody.extend(rows)
        tgroup += tbody

        return table

    def model_section(self, species, model):
        mid = f"sec_catalog_{species.id}_models_{model.id}"
        target = self.get_target(mid)
        section = nodes.section(ids=[mid])
        section += nodes.title(text=model.name)
        section += nodes.paragraph(text=model.description)
        section += nodes.rubric(text="Details")
        section += self.model_table(model)
        section += nodes.rubric(text="Populations")
        section += self.population_list(model)
        section += nodes.rubric(text="Citations")
        section += self.citation_list(model)
        section += nodes.rubric(text="Model parameters")
        section += self.model_parameter_table(species, model)
        return [target, section]

    def run(self):
        species = stdpopsim.get_species(self.arguments[0])
        sid = f"sec_catalog_{species.id}"
        species_target = self.get_target(sid)
        section = nodes.section(ids=[sid], names=[sid])
        section += nodes.title(text=species.name)
        section += self.species_summary_table(species)

        genome_section = nodes.section(ids=[f"sec_catalog_{species.id}_genome"])
        genome_section += nodes.title(text="Genome")
        genome_section += self.chromosomes_table(species)
        section += genome_section
        section += nodes.transition()

        maps_section = nodes.section(ids=[f"sec_catalog_{species.id}_genetic_maps"])
        maps_section += nodes.title(text="Genetic Maps")
        maps_section += self.genetic_maps_table(species)
        for gmap in species.genetic_maps:
            maps_section += self.genetic_map_section(species, gmap)
        section += maps_section
        section += nodes.transition()

        models_section = nodes.section(ids=[f"sec_catalog_{species.id}_models"])
        models_section += nodes.title(text="Models")
        # TODO we need some styling here to break this up visually. How does
        # sphinx autoclass do this? We should probably have a table of all the
        # models, and try to use the Sphinx styling to make it all work visually.
        # TODO add a table summarising the models with links to the detailed
        # descriptions
        for model in species.models:
            models_section += self.model_section(species, model)
        section += models_section

        return [species_target, section]


def setup(app):
    app.add_directive("speciescatalog", SpeciesCatalogDirective)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
