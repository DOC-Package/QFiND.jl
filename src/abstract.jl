abstract type DiscretizationMethod <: Function end
struct bsdo_discr <: DiscretizationMethod end
struct id_discr <: DiscretizationMethod end