"""
A sphinx extension that provides the "speciescatalog" directive. See the
Sphinx documentation for details on developing such extensions:

http://www.sphinx-doc.org/en/master/extdev/index.html#dev-extensions

This is based on the todo tutorial:
http://www.sphinx-doc.org/en/master/development/tutorials/todo.html
"""
import csv
import pathlib
import re

from docutils import nodes
from sphinx.util.docutils import SphinxDirective
from sphinx.util import logging

import stdpopsim
import os

logger = logging.getLogger(__name__)


class SpeciesCatalogDirective(SphinxDirective):

    required_arguments = 1

    def get_demographic_model_id(self, species, model):
        # It seems that we have to lowercase any IDs for sphinx. Need to
        # look into this, as it means we need to enforce IDs unique
        # without regard to case.
        return f"sec_catalog_{species.id}_models_{model.id}".lower()

    def get_genetic_map_id(self, species, genetic_map):
        return f"sec_catalog_{species.id}_genetic_maps_{genetic_map.id}".lower()

    def get_annotation_id(self, species, annotation):
        return f"sec_catalog_{species.id}_annotations_{annotation.id}".lower()

    def get_dfe_id(self, species, dfe):
        return f"sec_catalog_{species.id}_dfes_{dfe.id}".lower()

    def get_target(self, tid):
        """
        Returns a target node for the specified ID.

        This took a lot of digging through the Docutils and Sphinx sources to
        figure out: it appears that including the 'names' argument to the
        target node and calling note_explicit_target are essential to allow
        us to cross reference these targets from the rest of the documentation.
        """
        target = nodes.target("", "", ids=[tid], names=[tid])
        self.state.document.note_explicit_target(target)
        return target

    def make_field_list(self, data):

        field_list = nodes.field_list()
        for name, text, citations in data:
            field = nodes.field()
            field_name = nodes.field_name(text=name)
            field_body = nodes.field_body()
            para = nodes.paragraph(text=text)

            if citations is not None and len(citations) > 0:
                para += nodes.Text(" (")
                for i, citation in enumerate(citations):
                    text = f"{citation.author}, {citation.year}"
                    para += nodes.reference(
                        internal=False, refuri=citation.doi, text=text
                    )
                    if i != len(citations) - 1:
                        para += nodes.Text("; ")
                para += nodes.Text(")")

            field_body += para
            field += field_name
            field += field_body
            field_list += field

        return field_list

    def species_summary(self, species):
        data = [
            ("ID", species.id, None),
            ("Name", species.name, None),
            ("Common name", species.common_name, None),
            (
                "Generation time",
                species.generation_time,
                [
                    citation
                    for citation in species.citations
                    if stdpopsim.CiteReason.GEN_TIME in citation.reasons
                ],
            ),
            ("Ploidy", species.ploidy, None),
            (
                "Population size",
                species.population_size,
                [
                    citation
                    for citation in species.citations
                    if stdpopsim.CiteReason.POP_SIZE in citation.reasons
                ],
            ),
        ]
        return self.make_field_list(data)

    def model_summary(self, model):
        data = [
            ("ID", model.id, None),
            ("Description", model.description, None),
            ("Num populations", model.num_populations, None),
        ]
        return self.make_field_list(data)

    def citation_list(self, citable):
        bullet_list = nodes.bullet_list()
        for citation in citable.citations:
            list_item = nodes.list_item()

            para = nodes.paragraph(text=f"{citation.author}, {citation.year}. ")
            para += nodes.reference(
                internal=False, refuri=citation.doi, text=citation.doi
            )
            list_item += para
            bullet_list += list_item
        return bullet_list

    def model_parameter_table(self, species, model):
        path = pathlib.Path(f"parameter_tables/{species.id}/{model.id}.csv")
        if not path.exists():
            message = f"Skipping model parameters for {model.id} due to missing table"
            print(message)
            logger.warn(message)
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

    def population_table(self, model):
        table = nodes.table()
        tgroup = nodes.tgroup(cols=4)
        for _ in range(4):
            colspec = nodes.colspec(colwidth=1)
            tgroup.append(colspec)

        table += tgroup

        thead = nodes.thead()
        tgroup += thead
        row = nodes.row()
        entry = nodes.entry()
        entry += nodes.paragraph(text="Index")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text="ID")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text="Sampling time")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text="Description")
        row += entry
        thead.append(row)

        rows = []
        for index, population in enumerate(model.populations):
            row = nodes.row()
            rows.append(row)

            entry = nodes.entry()
            entry += nodes.paragraph(text=str(index))
            row += entry

            entry = nodes.entry()
            entry += nodes.paragraph(text=population.name)
            row += entry

            entry = nodes.entry()
            entry += nodes.paragraph(text=f"{population.default_sampling_time}")
            row += entry

            entry = nodes.entry()
            entry += nodes.paragraph(text=population.description)
            row += entry
        tbody = nodes.tbody()
        tbody.extend(rows)
        tgroup += tbody
        return table

    def genetic_map_section(self, species, genetic_map):
        map_id = self.get_genetic_map_id(species, genetic_map)
        target = self.get_target(map_id)
        section = nodes.section(ids=[map_id])
        section += nodes.title(text=genetic_map.id)
        section += nodes.paragraph(text=genetic_map.long_description)
        section += nodes.rubric(text="Citations")
        section += self.citation_list(genetic_map)
        return [target, section]

    def genetic_maps_table(self, species):
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
        entry += nodes.paragraph(text="ID")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text="Year")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text="Description")
        row += entry

        thead.append(row)

        rows = []
        for genetic_map in species.genetic_maps:
            row = nodes.row()
            rows.append(row)

            map_id = self.get_genetic_map_id(species, genetic_map)
            entry = nodes.entry()
            para = nodes.paragraph()
            entry += para
            para += nodes.reference(internal=True, refid=map_id, text=genetic_map.id)
            row += entry

            entry = nodes.entry()
            entry += nodes.paragraph(text=genetic_map.citations[0].year)
            row += entry

            entry = nodes.entry()
            para = nodes.paragraph()
            entry += nodes.paragraph(text=genetic_map.description)
            row += entry

        tbody = nodes.tbody()
        tbody.extend(rows)
        tgroup += tbody

        return table

    def chromosomes_table(self, species):

        table = nodes.table()
        tgroup = nodes.tgroup(cols=5)
        for _ in range(5):
            colspec = nodes.colspec(colwidth=1)
            tgroup.append(colspec)
        table += tgroup

        thead = nodes.thead()
        tgroup += thead
        row = nodes.row()
        entry = nodes.entry()
        entry += nodes.paragraph(text="ID")
        row += entry

        entry = nodes.entry()
        entry += nodes.paragraph(text="Ploidy")
        row += entry

        entry = nodes.entry()
        entry += nodes.paragraph(text="Length")
        row += entry

        entry = nodes.entry()
        entry += nodes.paragraph(text="Recombination rate")
        row += entry

        entry = nodes.entry()
        entry += nodes.paragraph(text="Mutation rate")
        row += entry

        thead.append(row)

        rows = []
        for chrom in species.genome.chromosomes:
            row = nodes.row()
            entry = nodes.entry()
            entry += nodes.paragraph(text=chrom.id)
            row += entry

            entry = nodes.entry()
            entry += nodes.paragraph(text="{:d}".format(chrom.ploidy))
            row += entry

            entry = nodes.entry()
            entry += nodes.paragraph(text="{:d}".format(chrom.length))
            row += entry

            entry = nodes.entry()
            entry += nodes.paragraph(text="{:g}".format(chrom.recombination_rate))
            row += entry

            entry = nodes.entry()
            entry += nodes.paragraph(text="{:g}".format(chrom.mutation_rate))
            row += entry

            rows.append(row)
        tbody = nodes.tbody()
        tbody.extend(rows)
        tgroup += tbody

        return table

    def models_table(self, species):
        table = nodes.table()
        tgroup = nodes.tgroup(cols=2)
        for _ in range(2):
            colspec = nodes.colspec(colwidth=1)
            tgroup.append(colspec)
        table += tgroup

        thead = nodes.thead()
        tgroup += thead
        row = nodes.row()
        entry = nodes.entry()
        entry += nodes.paragraph(text="ID")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text="Description")
        row += entry

        thead.append(row)

        rows = []
        for model in species.demographic_models:
            row = nodes.row()
            rows.append(row)

            mid = self.get_demographic_model_id(species, model)
            entry = nodes.entry()
            para = nodes.paragraph()
            entry += para
            para += nodes.reference(internal=True, refid=mid, text=model.id)
            row += entry

            entry = nodes.entry()
            entry += nodes.paragraph(text=model.description)
            row += entry

        tbody = nodes.tbody()
        tbody.extend(rows)
        tgroup += tbody

        return table

    def model_section(self, species, model):
        mid = self.get_demographic_model_id(species, model)
        target = self.get_target(mid)
        section = nodes.section(ids=[mid])
        section += nodes.title(text=model.description)
        section += nodes.paragraph(text=model.long_description)
        section += nodes.rubric(text="Details")
        section += self.model_summary(model)
        section += nodes.rubric(text="Populations")
        section += self.population_table(model)
        section += nodes.rubric(text="Citations")
        section += self.citation_list(model)
        section += nodes.rubric(text="Demographic Model parameters")
        section += self.model_parameter_table(species, model)
        return [target, section]

    def model_image(self, species, model):
        import demesdraw
        import matplotlib.pyplot as plt

        mid = self.get_demographic_model_id(species, model)
        _, ax = plt.subplots(1, 1, figsize=(4.5, 4), tight_layout=True)
        # Conversion into demes object for easier plotting
        graph = model.model.to_demes()
        ax = demesdraw.tubes(graph, ax=ax, log_time=True)
        ax.set_title(f"{model.id}", fontsize=10)
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=45,
            ha="right",
            rotation_mode="anchor",
            fontsize=10,
        )
        ax.set_ylabel("Time (generations)", fontsize=10)
        os.makedirs(f"parameter_images/{species.id}/", exist_ok=True)
        img_name = f"parameter_images/{species.id}/{mid}.png"
        plt.savefig(img_name, dpi=150)
        plt.close()
        section = nodes.image(uri=img_name)
        return section

    def annotation_section(self, species, annotation):
        map_id = self.get_annotation_id(species, annotation)
        target = self.get_target(map_id)
        section = nodes.section(ids=[map_id])
        section += nodes.title(text=annotation.id)
        section += nodes.paragraph(text=annotation.description)
        section += nodes.rubric(text="Citations")
        section += self.citation_list(annotation)
        return [target, section]

    def annotation_table(self, species):
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
        entry += nodes.paragraph(text="ID")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text="Year")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text="Description")
        row += entry

        thead.append(row)

        rows = []
        for an in species.annotations:
            row = nodes.row()
            rows.append(row)

            map_id = self.get_annotation_id(species, an)
            entry = nodes.entry()
            para = nodes.paragraph()
            entry += para
            para += nodes.reference(internal=True, refid=map_id, text=an.id)
            row += entry

            entry = nodes.entry()
            entry += nodes.paragraph(text=an.citations[0].year)
            row += entry

            entry = nodes.entry()
            para = nodes.paragraph()
            entry += nodes.paragraph(text=an.description)
            row += entry

        tbody = nodes.tbody()
        tbody.extend(rows)
        tgroup += tbody

        return table

    def dfe_section(self, species, dfe):
        dfe_id = self.get_dfe_id(species, dfe)
        target = self.get_target(dfe_id)
        section = nodes.section(ids=[dfe_id])
        section += nodes.title(text=dfe.id)
        section += nodes.paragraph(text=dfe.description)
        section += nodes.rubric(text="Citations")
        section += self.citation_list(dfe)
        section += nodes.rubric(text="DFE parameters")
        section += self.dfe_parameter_table(species, dfe)
        return [target, section]

    def dfe_parameter_table(self, species, dfe):
        table = nodes.table()
        # proportion, distribution type, arguments, dominance
        tgroup = nodes.tgroup(cols=4)
        for _ in range(4):
            colspec = nodes.colspec(colwidth=1)
            tgroup.append(colspec)

        table += tgroup

        thead = nodes.thead()
        tgroup += thead
        row = nodes.row()
        entry = nodes.entry()
        entry += nodes.paragraph(text="Proportion of mutations")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text="Distribution type")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text="Parameters")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text="Dominance")
        row += entry
        thead.append(row)

        def parse_mutation_type(mt):
            dt = mt.distribution_type
            da = mt.distribution_args
            if dt == "s":
                # currently we only use this for rlnorm and unif
                bits = re.split("[ ()*]", da[0])
                assert bits[0] == "return", "Unrecognized mutation type string"
                if "rlnorm" in bits:
                    dt = "ln" if bits[1] == "-1" else "lp"
                    da = (float(bits[bits.index("rlnorm") + 2]), float(bits[-2]))
                elif "runif" in bits:
                    dt = "u"
                    da = (float(bits[6]), float(bits[10]))
                else:
                    raise ValueError("Unrecognized mutation type string")

            if dt == "f":
                label = "Fixed s"
                params = f"s = {da[0]:.3f}"
            elif dt == "e":
                label = "Exponential"
                params = f"mean = {da[0]:.3f}"
            elif dt == "g":
                label = "Gamma"
                params = f"mean = {da[0]:.3f}, " f"shape = {da[1]:.3f}"
            elif dt == "n":
                label = "Normal"
                params = f"mean = {da[0]:.3f}, " f"sd = {da[1]:.3f}"
            elif dt == "w":
                label = "Weibull"
                params = f"scale = {da[0]:.3f}, " f"shape = {da[1]:.3f}"
            elif dt == "u":
                label = "Uniform"
                params = f"min = {da[0]:.3f}, " f"max = {da[1]:.3f}"
            elif dt == "lp":
                label = "Positive LogNormal"
                params = f"meanlog = {da[0]:.3f}, " f"sdlog = {da[1]:.3f}"
            elif dt == "ln":
                label = "Negative LogNormal"
                params = f"meanlog = {da[0]:.3f}, " f"sdlog = {da[1]:.3f}"

            if mt.dominance_coeff is None:
                h_list = []
                for j in range(len(mt.dominance_coeff_list)):
                    slist = (
                        ["-Inf"]
                        + [f"{s:.3f}" for s in mt.dominance_coeff_breaks]
                        + ["Inf"]
                    )
                    if dt in ["e", "g", "w"]:
                        if da[0] < 0:
                            slist[0] = "0"
                        else:
                            slist[-1] = "0"
                    elif dt == "lp":
                        slist[0] = "0"
                    elif dt == "ln":
                        slist[-1] = "0"
                    mins = slist[j]
                    maxs = slist[j + 1]
                    h_list.append(
                        f"{mt.dominance_coeff_list[j]:.3f} on " f"[{mins}, " f"{maxs})"
                    )
                h = "h = " + "; ".join(h_list)
            else:
                h = f"h = {mt.dominance_coeff:.3f}"
            return label, params, h

        rows = []
        for mt, prop in zip(dfe.mutation_types, dfe.proportions):
            if prop > 0:
                row = nodes.row()
                rows.append(row)
                entry = nodes.entry()
                entry += nodes.paragraph(text=f"{100*prop:.1f}%")
                row += entry
                for label in parse_mutation_type(mt):
                    entry = nodes.entry()
                    entry += nodes.paragraph(text=label)
                    row += entry

        tbody = nodes.tbody()
        tbody.extend(rows)
        tgroup += tbody

        return table

    def dfes_table(self, species):
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
        entry += nodes.paragraph(text="ID")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text="Year")
        row += entry
        entry = nodes.entry()
        entry += nodes.paragraph(text="Description")
        row += entry

        thead.append(row)

        rows = []

        for dfe in species.dfes:
            row = nodes.row()
            rows.append(row)

            dfe_id = self.get_dfe_id(species, dfe)
            entry = nodes.entry()
            para = nodes.paragraph()
            entry += para
            para += nodes.reference(internal=True, refid=dfe_id, text=dfe.id)
            row += entry

            entry = nodes.entry()
            entry += nodes.paragraph(text=dfe.citations[0].year)
            row += entry

            entry = nodes.entry()
            para = nodes.paragraph()
            entry += nodes.paragraph(text=dfe.description)
            row += entry

        tbody = nodes.tbody()
        tbody.extend(rows)
        tgroup += tbody

        return table

    def run(self):
        # species:
        species = stdpopsim.get_species(self.arguments[0])
        sid = f"sec_catalog_{species.id}"
        species_target = self.get_target(sid.lower())
        section = nodes.section(ids=[sid], names=[sid])
        section += nodes.title(text=species.name)
        section += self.species_summary(species)
        # genomes:
        genome_section = nodes.section(ids=[f"sec_catalog_{species.id}_genome"])
        genome_section += nodes.title(text="Genome")
        if species.genome.bacterial_recombination is True:
            genome_section += self.make_field_list(
                [
                    (
                        "Bacterial recombination with tract length range",
                        species.genome.range_gene_conversion_lengths,
                        None,
                    )
                ]
            )
        elif species.genome.mean_gene_conversion_fraction > 0:
            genome_section += self.make_field_list(
                [
                    (
                        "Mean gene conversion fraction",
                        species.genome.mean_gene_conversion_fraction,
                        None,
                    ),
                    (
                        "Range gene conversion lengths",
                        species.genome.range_gene_conversion_lengths,
                        None,
                    ),
                ]
            )
        genome_section += self.chromosomes_table(species)
        genome_section += nodes.paragraph(
            text="Mutation and recombination rates "
            "are in units of per bp and per generation."
        )
        section += genome_section
        section += nodes.transition()
        # genetic maps:
        if len(species.genetic_maps) > 0:
            map_id = f"sec_catalog_{species.id}_genetic_maps"
            maps_section = nodes.section(ids=[map_id])
            maps_section += nodes.title(text="Genetic Maps")
            maps_section += self.genetic_maps_table(species)
            for gmap in species.genetic_maps:
                maps_section += self.genetic_map_section(species, gmap)
            section += self.get_target(map_id.lower())
            section += maps_section
            section += nodes.transition()
        # demographic models:
        if len(species.demographic_models) > 0:
            map_id = f"sec_catalog_{species.id}_models"
            models_section = nodes.section(ids=[map_id])
            models_section += nodes.title(text="Demographic Models")
            models_section += self.models_table(species)
            for i, model in enumerate(species.demographic_models):
                models_section += self.model_section(species, model)
                models_section += self.model_image(species, model)
                if i < len(species.demographic_models) - 1:
                    models_section += nodes.transition()
            section += self.get_target(map_id.lower())
            section += models_section
        # annotation:
        if len(species.annotations) > 0:
            map_id = f"sec_catalog_{species.id}_annotations"
            annot_section = nodes.section(ids=[map_id])
            annot_section += nodes.title(text="Annotations")
            annot_section += self.annotation_table(species)
            for an in species.annotations:
                annot_section += self.annotation_section(species, an)
            section += self.get_target(map_id.lower())
            section += annot_section
            section += nodes.transition()
        # DFE:
        if len(species.dfes) > 0:
            map_id = f"sec_catalog_{species.id}_dfe"
            dfes_section = nodes.section(ids=[map_id])
            dfes_section += nodes.title(text="Distribution of Fitness Effects (DFEs)")
            dfes_section += self.dfes_table(species)
            for i, dfe in enumerate(species.dfes):
                dfes_section += self.dfe_section(species, dfe)
            section += self.get_target(map_id.lower())
            section += dfes_section
            section += nodes.transition()
        return [species_target, section]


def setup(app):
    app.add_directive("speciescatalog", SpeciesCatalogDirective)

    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
