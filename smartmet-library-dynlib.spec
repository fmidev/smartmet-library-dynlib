%bcond_with tests
%define DIRNAME dynlib
%define LIBNAME smartmet-%{DIRNAME}
%define SPECNAME smartmet-library-%{DIRNAME}

Summary: dynamic-meteorology feature detection for SmartMet
Name: %{SPECNAME}
Version: 26.4.17
Release: 1%{?dist}.fmi
License: MIT
Group: Development/Libraries
URL: https://github.com/fmidev/smartmet-library-dynlib
Source0: %{name}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)

BuildRequires: gcc-c++
BuildRequires: gcc-gfortran
BuildRequires: lapack-devel
BuildRequires: make
BuildRequires: rpm-build
BuildRequires: smartmet-library-macgyver-devel
BuildRequires: smartmet-utils-devel >= 26.2.4

Requires: blas
Requires: lapack
Requires: libgfortran

Provides: %{SPECNAME}

%description
Feature-detection routines for dynamic meteorology, backed by the
dynlib Fortran library of Clemens Spensberger (University of Bergen).
Exposes front, jet, trough, cyclone, Rossby wave breaking, blocking,
and atmospheric river detection to C++ callers in the SmartMet
ecosystem. The upstream Fortran kernels are vendored under
third_party/dynlib/ with their original MIT licence preserved.

%prep
rm -rf $RPM_BUILD_ROOT
rm -rf %{SPECNAME}

%setup -q -n %{SPECNAME}

%build
make %{_smp_mflags}
%if %{with tests}
make test %{_smp_mflags}
%endif

%install
%makeinstall

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig
%postun -p /sbin/ldconfig

%files
%defattr(0775,root,root,0775)
%{_libdir}/libsmartmet-%{DIRNAME}.so

%package -n %{SPECNAME}-devel
Summary: development files for smartmet-library-dynlib
Provides: %{SPECNAME}-devel
Requires: %{SPECNAME} = %{version}-%{release}
Requires: lapack-devel
Requires: smartmet-library-macgyver-devel

%description -n %{SPECNAME}-devel
Headers and development files for smartmet-library-dynlib.

%files -n %{SPECNAME}-devel
%defattr(0664,root,root,0775)
%{_includedir}/smartmet/%{DIRNAME}

%changelog
* Fri Apr 17 2026 Mika Heiskanen <mika.heiskanen@fmi.fi> 26.4.17-1.fmi
- Initial packaging: vendored dynlib Fortran sources plus ISO_C_BINDING
  shim, exposing Jenkner-style maximum-gradient front detection to
  C++ callers
