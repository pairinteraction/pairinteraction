{{ name | escape | underline }}
Full {{ objtype }} path: {{ fullname | escape }}

.. workaround to ignore generic classes like AtomSQDT[int]
{% set filtered_classes = [] %}
{% for item in classes %}
   {% if not item.endswith(']') %}
      {% set _ = filtered_classes.append(item) %}
   {% endif %}
{%- endfor %}
{% set classes = filtered_classes %}

.. currentmodule:: {{ fullname }}

{% if attributes %}
.. rubric:: {{ _('Module Attributes') }}

.. autosummary::
   {% for item in attributes %}
      {{ item }}
   {%- endfor %}
{% endif %}

{% if functions %}
.. rubric:: {{ _('Functions') }}

.. autosummary::
   :toctree: .
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
{% endif %}

{% if classes %}
.. rubric:: {{ _('Classes') }}

.. .. inheritance-diagram:: {{ fullname }}
..     :parts: 1

.. autosummary::
   :toctree: .
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
{% endif %}

{% if exceptions %}
.. rubric:: {{ _('Exceptions') }}

.. autosummary::
   :toctree: .
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
{% endif %}

{% if modules %}
.. rubric:: Modules

.. autosummary::
   :toctree: .
   :recursive:
   {% for item in modules %}
      {{ item }}
   {%- endfor %}
{% endif %}

.. rubric:: Module description

.. automodule:: {{ fullname }}
