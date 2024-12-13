{{ name | escape | underline }}
Full {{ objtype }} path: {{ fullname | escape }}

.. inheritance-diagram:: {{ fullname }}
    :parts: 1

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
   :undoc-members:
   :class-doc-from: init
