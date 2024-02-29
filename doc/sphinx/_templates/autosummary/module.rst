{{ name | escape | underline }}

{% set pydantic_models = [] %}
{% for item in members %}
   {% if item.startswith('BaseModel') or item.startswith('Model') or item.startswith('BaseParameter') or item.startswith('Parameter') %}
      {% set _ = pydantic_models.append(item) %}
   {% endif %}
{%- endfor %}

- full name: {{ fullname | escape }}
- type: {{ objtype }}

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

.. inheritance-diagram:: {{ fullname }}
    :parts: 1

.. autosummary::
   :toctree: .
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
{% endif %}

{% if pydantic_models %}
.. rubric:: {{ _('Pydantic Models') }}

.. inheritance-diagram:: {{ fullname }}
    :parts: 1

.. autosummary::
   :toctree: .
   {% for item in pydantic_models %}
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
   :show-inheritance:
