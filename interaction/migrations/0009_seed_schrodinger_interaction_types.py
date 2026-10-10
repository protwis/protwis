"""Seed every interaction type the Schrodinger imports write."""
# build_structures does not compute ligand interactions, so nothing else in a
# build creates these types. These rows are exactly
# schrodinger_import.required_slugs(): the targets of interaction_type_map.yaml
# outside the excluded families, plus polar_backbone, the backbone override.
#
# Eleven of them are types the legacy calculation also computes, and keep the
# slug, name, type and direction it gives them, so a database built either way
# names them alike. Three are new:
#
# * halogen_protein and metal_coordination_protein. halogen_protein is named
#   "halogen contact", not "halogen bond": the producer's criterion has loose
#   angle limits, so the rows are contacts rather than proven halogen bonds.
# * covalent, one row per ligand-receptor bond; a row says a bond exists, not
#   its order. Its type is ``covalent`` and not ``hidden``: pages leave out
#   hidden types, and the rows are imported so that the pages show them.
#
# Existing rows are never modified, except that an existing halogen_protein row
# is renamed to "halogen contact". The reverse operation is a deliberate no-op:
# deleting interaction types would cascade to every ResidueFragmentInteraction
# that uses them.

from django.db import migrations


# The order is the order rows are created in, so it sets the ids new rows get
# (the ligand page lists a ligand's interaction rows in type id order).
SEEDED_TYPES = (
    # slug, name, type, direction
    ("aro_ion_protein", "aromatic (pi-cation)", "aromatic", "protein"),
    ("halogen_protein", "halogen contact", "polar", ""),
    ("metal_coordination_protein", "metal coordination", "polar", ""),
    ("covalent", "covalent bond", "covalent", ""),
    ("acc", "accessible", "hidden", ""),
    ("aro_ef_protein", "aromatic (edge-to-face)", "aromatic", "protein"),
    ("aro_ff", "aromatic (face-to-face)", "aromatic", "none"),
    ("hyd", "hydrophobic", "hydrophobic", ""),
    ("polar_acceptor_protein", "polar (hydrogen bond)", "polar", "protein"),
    ("polar_backbone", "polar (hydrogen bond with backbone)", "polar", "protein"),
    ("polar_donor_protein", "polar (hydrogen bond)", "polar", "protein"),
    ("polar_double_neg_protein", "polar (charge-charge)", "polar", ""),
    ("polar_double_pos_protein", "polar (charge-charge)", "polar", ""),
    ("Van der Waals", "Van der Waals", "waals", ""),
)


def seed(apps, schema_editor):
    InteractionType = apps.get_model("interaction", "ResidueFragmentInteractionType")
    for slug, name, type_, direction in SEEDED_TYPES:
        row, created = InteractionType.objects.get_or_create(
            slug=slug, defaults={"name": name, "type": type_, "direction": direction}
        )
        # An existing database may name this slug "halogen bond"; the name is
        # the only field corrected.
        if not created and slug == "halogen_protein" and row.name != name:
            row.name = name
            row.save(update_fields=["name"])


class Migration(migrations.Migration):
    dependencies = [
        ("interaction", "0008_auto_20260921_1803"),
    ]

    operations = [
        migrations.RunPython(seed, migrations.RunPython.noop),
    ]
