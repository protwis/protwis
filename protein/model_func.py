from protein.models import ProteinFamily

def get_current_classless_parent_gpcr_slugs(old_slug_dict, new_slug_set):
    """ This is a solution for new builds with change of the definition of CLASSLESS_PARENT_GPCR_SLUGS

        old_slug_dict:  a dict with the deprecated classless parent GPCR family slugs. 
                        Format: {CLASSLESS_PARENT_GPCR_SLUG:CLASSLESS_PARENT_GPCR_NAME}

        new_slug_set:   The new CLASSLESS_PARENT_GPCR_SLUGS constant value. A set that contains the new
                        classless parent GPCR family slugs.

        Example:    from protein.model_func import get_current_classless_parent_gpcr_slugs
                    CLASSLESS_PARENT_GPCR_SLUGS = get_current_classless_parent_gpcr_slugs({'008':'Other GPCRs'}, {'010'})
    
        If slugs and family names in old_slug_dict do not match with the ones in the database or they are 
        missing in the database, returns new_slug_set. It also returns new_slug_set if a slug of new_slug_set 
        not present in old_slug_dict is present in the database. Otherwise, returns a set that contains
        the keys of old_slug_dict.

        The function must be imported and run in protein.models after class ProteinFamily is defined.

    """
    old_slug_dict_keys = old_slug_dict.keys()
    old_slug_dict_keys_set = set(old_slug_dict_keys)

    # Make a set of the old slugs that are currently present in the database and check if their names match
    # with the old ones
    q = ProteinFamily.objects.filter(slug__in=list(old_slug_dict_keys)).values('slug','name')
    current_db_old_slugs_set = set()
    for pf in q:
        current_db_old_slugs_set.add(pf['slug'])
        if pf['name'] != old_slug_dict[pf['slug']]:
            return new_slug_set

    # Check if any new slugs are present in the database
    new_slugs = new_slug_set - old_slug_dict_keys_set
    if len(ProteinFamily.objects.filter(slug__in=list(new_slugs)).values_list('slug', flat=True)) > 0:
        return new_slug_set

    # Check if any old slugs are missing in the database
    if current_db_old_slugs_set != old_slug_dict_keys_set:
        return new_slug_set

    return old_slug_dict_keys_set

