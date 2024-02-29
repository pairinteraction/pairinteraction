{{ name | escape | underline }}

.. inheritance-diagram:: {{ fullname }}
    :parts: 1

- full name: {{ fullname | escape }}
- type: {{ objtype }}

.. currentmodule:: {{ module }}

{% if methods %}
.. rubric:: Class Methods

.. autosummary::
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
{% endif %}

{% if attributes %}
.. rubric:: Class Attributes and Properties

.. autosummary::
{% for item in attributes %}
   ~{{ name }}.{{ item }}
{%- endfor %}
{% endif %}


.. autoclass:: {{ fullname }}
   :members:
   :member-order: bysource
   :inherited-members:
   :show-inheritance:
