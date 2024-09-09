{{ name | escape | underline }}

.. currentmodule:: {{ module }}

{% if objtype.startswith('pydantic_model') %}
.. inheritance-diagram:: {{ fullname }}
    :parts: 1
{% endif %}

- full name: {{ fullname | escape }}
- type: {{ objtype }}

.. auto{{ objtype }}:: {{ objname }}
