%define skip_python2 1
Name:           python-sphinxcontrib-mermaid
Version:        0.7.1
Release:        0
Summary:        Mermaid diagrams in yours sphinx powered docs
License:        BSD-2-Clause
Group:          Development/Languages/Python
URL:            https://github.com/mgaitan/sphinxcontrib-mermaid
Source:         https://files.pythonhosted.org/packages/source/s/sphinxcontrib-mermaid/sphinxcontrib-mermaid-%{version}.tar.gz
BuildRequires:  python3-setuptools
BuildRequires:  fdupes
BuildRequires:  python-rpm-macros
Requires:       python-Sphinx
BuildArch:      noarch
%python_subpackages

%description
This extension allows you to embed Mermaid graphs in your documents,
including general flowcharts, sequence diagrams, gantt diagrams and more.

%prep
%setup -q -n sphinxcontrib-mermaid-%{version}

%build
%python_build

%install
%python_install
%python_expand %fdupes %{buildroot}%{$python_sitelib}

%files %{python_files}
%doc README.rst CHANGELOG.rst
%license LICENSE.rst
%{python_sitelib}/*

%changelog
