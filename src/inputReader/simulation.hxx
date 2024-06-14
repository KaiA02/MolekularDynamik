// Copyright (c) 2005-2014 Code Synthesis Tools CC
//
// This program was generated by CodeSynthesis XSD, an XML Schema to
// C++ data binding compiler.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 2 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
//
// In addition, as a special exception, Code Synthesis Tools CC gives
// permission to link this program with the Xerces-C++ library (or with
// modified versions of Xerces-C++ that use the same license as Xerces-C++),
// and distribute linked combinations including the two. You must obey
// the GNU General Public License version 2 in all respects for all of
// the code used other than Xerces-C++. If you modify this copy of the
// program, you may extend this exception to your version of the program,
// but you are not obligated to do so. If you do not wish to do so, delete
// this exception statement from your version.
//
// Furthermore, Code Synthesis Tools CC makes a special exception for
// the Free/Libre and Open Source Software (FLOSS) which is described
// in the accompanying FLOSSE file.
//

#ifndef SIMULATION_HXX
#define SIMULATION_HXX

#ifndef XSD_CXX11
#define XSD_CXX11
#endif

#ifndef XSD_USE_CHAR
#define XSD_USE_CHAR
#endif

#ifndef XSD_CXX_TREE_USE_CHAR
#define XSD_CXX_TREE_USE_CHAR
#endif

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/config.hxx>

#if (XSD_INT_VERSION != 4000000L)
#error XSD runtime version mismatch
#endif

#include <xsd/cxx/pre.hxx>

#include <xsd/cxx/xml/char-utf8.hxx>

#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/elements.hxx>
#include <xsd/cxx/tree/types.hxx>

#include <xsd/cxx/xml/error-handler.hxx>

#include <xsd/cxx/xml/dom/auto-ptr.hxx>

#include <xsd/cxx/tree/parsing.hxx>
#include <xsd/cxx/tree/parsing/byte.hxx>
#include <xsd/cxx/tree/parsing/unsigned-byte.hxx>
#include <xsd/cxx/tree/parsing/short.hxx>
#include <xsd/cxx/tree/parsing/unsigned-short.hxx>
#include <xsd/cxx/tree/parsing/int.hxx>
#include <xsd/cxx/tree/parsing/unsigned-int.hxx>
#include <xsd/cxx/tree/parsing/long.hxx>
#include <xsd/cxx/tree/parsing/unsigned-long.hxx>
#include <xsd/cxx/tree/parsing/boolean.hxx>
#include <xsd/cxx/tree/parsing/float.hxx>
#include <xsd/cxx/tree/parsing/double.hxx>
#include <xsd/cxx/tree/parsing/decimal.hxx>

namespace xml_schema
{
  // anyType and anySimpleType.
  //
  typedef ::xsd::cxx::tree::type type;
  typedef ::xsd::cxx::tree::simple_type< char, type > simple_type;
  typedef ::xsd::cxx::tree::type container;

  // 8-bit
  //
  typedef signed char byte;
  typedef unsigned char unsigned_byte;

  // 16-bit
  //
  typedef short short_;
  typedef unsigned short unsigned_short;

  // 32-bit
  //
  typedef int int_;
  typedef unsigned int unsigned_int;

  // 64-bit
  //
  typedef long long long_;
  typedef unsigned long long unsigned_long;

  // Supposed to be arbitrary-length integral types.
  //
  typedef long long integer;
  typedef long long non_positive_integer;
  typedef unsigned long long non_negative_integer;
  typedef unsigned long long positive_integer;
  typedef long long negative_integer;

  // Boolean.
  //
  typedef bool boolean;

  // Floating-point types.
  //
  typedef float float_;
  typedef double double_;
  typedef double decimal;

  // String types.
  //
  typedef ::xsd::cxx::tree::string< char, simple_type > string;
  typedef ::xsd::cxx::tree::normalized_string< char, string > normalized_string;
  typedef ::xsd::cxx::tree::token< char, normalized_string > token;
  typedef ::xsd::cxx::tree::name< char, token > name;
  typedef ::xsd::cxx::tree::nmtoken< char, token > nmtoken;
  typedef ::xsd::cxx::tree::nmtokens< char, simple_type, nmtoken > nmtokens;
  typedef ::xsd::cxx::tree::ncname< char, name > ncname;
  typedef ::xsd::cxx::tree::language< char, token > language;

  // ID/IDREF.
  //
  typedef ::xsd::cxx::tree::id< char, ncname > id;
  typedef ::xsd::cxx::tree::idref< char, ncname, type > idref;
  typedef ::xsd::cxx::tree::idrefs< char, simple_type, idref > idrefs;

  // URI.
  //
  typedef ::xsd::cxx::tree::uri< char, simple_type > uri;

  // Qualified name.
  //
  typedef ::xsd::cxx::tree::qname< char, simple_type, uri, ncname > qname;

  // Binary.
  //
  typedef ::xsd::cxx::tree::buffer< char > buffer;
  typedef ::xsd::cxx::tree::base64_binary< char, simple_type > base64_binary;
  typedef ::xsd::cxx::tree::hex_binary< char, simple_type > hex_binary;

  // Date/time.
  //
  typedef ::xsd::cxx::tree::time_zone time_zone;
  typedef ::xsd::cxx::tree::date< char, simple_type > date;
  typedef ::xsd::cxx::tree::date_time< char, simple_type > date_time;
  typedef ::xsd::cxx::tree::duration< char, simple_type > duration;
  typedef ::xsd::cxx::tree::gday< char, simple_type > gday;
  typedef ::xsd::cxx::tree::gmonth< char, simple_type > gmonth;
  typedef ::xsd::cxx::tree::gmonth_day< char, simple_type > gmonth_day;
  typedef ::xsd::cxx::tree::gyear< char, simple_type > gyear;
  typedef ::xsd::cxx::tree::gyear_month< char, simple_type > gyear_month;
  typedef ::xsd::cxx::tree::time< char, simple_type > time;

  // Entity.
  //
  typedef ::xsd::cxx::tree::entity< char, ncname > entity;
  typedef ::xsd::cxx::tree::entities< char, simple_type, entity > entities;

  typedef ::xsd::cxx::tree::content_order content_order;
  // Flags and properties.
  //
  typedef ::xsd::cxx::tree::flags flags;
  typedef ::xsd::cxx::tree::properties< char > properties;

  // Parsing/serialization diagnostics.
  //
  typedef ::xsd::cxx::tree::severity severity;
  typedef ::xsd::cxx::tree::error< char > error;
  typedef ::xsd::cxx::tree::diagnostics< char > diagnostics;

  // Exceptions.
  //
  typedef ::xsd::cxx::tree::exception< char > exception;
  typedef ::xsd::cxx::tree::bounds< char > bounds;
  typedef ::xsd::cxx::tree::duplicate_id< char > duplicate_id;
  typedef ::xsd::cxx::tree::parsing< char > parsing;
  typedef ::xsd::cxx::tree::expected_element< char > expected_element;
  typedef ::xsd::cxx::tree::unexpected_element< char > unexpected_element;
  typedef ::xsd::cxx::tree::expected_attribute< char > expected_attribute;
  typedef ::xsd::cxx::tree::unexpected_enumerator< char > unexpected_enumerator;
  typedef ::xsd::cxx::tree::expected_text_content< char > expected_text_content;
  typedef ::xsd::cxx::tree::no_prefix_mapping< char > no_prefix_mapping;

  // Error handler callback interface.
  //
  typedef ::xsd::cxx::xml::error_handler< char > error_handler;

  // DOM interaction.
  //
  namespace dom
  {
    // Automatic pointer for DOMDocument.
    //
    using ::xsd::cxx::xml::dom::unique_ptr;

#ifndef XSD_CXX_TREE_TREE_NODE_KEY__XML_SCHEMA
#define XSD_CXX_TREE_TREE_NODE_KEY__XML_SCHEMA
    // DOM user data key for back pointers to tree nodes.
    //
    const XMLCh* const tree_node_key = ::xsd::cxx::tree::user_data_keys::node;
#endif
  }
}

// Forward declarations.
//
class simulation;
class input;
class output;
class config;
class particles;
class cuboids;
class disk;

#include <memory>    // ::std::unique_ptr
#include <limits>    // std::numeric_limits
#include <algorithm> // std::binary_search
#include <utility>   // std::move

#include <xsd/cxx/xml/char-utf8.hxx>

#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/elements.hxx>
#include <xsd/cxx/tree/containers.hxx>
#include <xsd/cxx/tree/list.hxx>

#include <xsd/cxx/xml/dom/parsing-header.hxx>

class simulation: public ::xml_schema::type
{
  public:
  // input
  //
  typedef ::input input_type;
  typedef ::xsd::cxx::tree::traits< input_type, char > input_traits;

  const input_type&
  input () const;

  input_type&
  input ();

  void
  input (const input_type& x);

  void
  input (::std::unique_ptr< input_type > p);

  // output
  //
  typedef ::output output_type;
  typedef ::xsd::cxx::tree::traits< output_type, char > output_traits;

  const output_type&
  output () const;

  output_type&
  output ();

  void
  output (const output_type& x);

  void
  output (::std::unique_ptr< output_type > p);

  // config
  //
  typedef ::config config_type;
  typedef ::xsd::cxx::tree::traits< config_type, char > config_traits;

  const config_type&
  config () const;

  config_type&
  config ();

  void
  config (const config_type& x);

  void
  config (::std::unique_ptr< config_type > p);

  // Constructors.
  //
  simulation (const input_type&,
              const output_type&,
              const config_type&);

  simulation (::std::unique_ptr< input_type >,
              ::std::unique_ptr< output_type >,
              ::std::unique_ptr< config_type >);

  simulation (const ::xercesc::DOMElement& e,
              ::xml_schema::flags f = 0,
              ::xml_schema::container* c = 0);

  simulation (const simulation& x,
              ::xml_schema::flags f = 0,
              ::xml_schema::container* c = 0);

  virtual simulation*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  simulation&
  operator= (const simulation& x);

  virtual 
  ~simulation ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< input_type > input_;
  ::xsd::cxx::tree::one< output_type > output_;
  ::xsd::cxx::tree::one< config_type > config_;
};

class input: public ::xml_schema::type
{
  public:
  // tStart
  //
  typedef ::xml_schema::double_ tStart_type;
  typedef ::xsd::cxx::tree::traits< tStart_type, char, ::xsd::cxx::tree::schema_type::double_ > tStart_traits;

  const tStart_type&
  tStart () const;

  tStart_type&
  tStart ();

  void
  tStart (const tStart_type& x);

  // tEnd
  //
  typedef ::xml_schema::double_ tEnd_type;
  typedef ::xsd::cxx::tree::traits< tEnd_type, char, ::xsd::cxx::tree::schema_type::double_ > tEnd_traits;

  const tEnd_type&
  tEnd () const;

  tEnd_type&
  tEnd ();

  void
  tEnd (const tEnd_type& x);

  // deltaT
  //
  typedef ::xml_schema::double_ deltaT_type;
  typedef ::xsd::cxx::tree::traits< deltaT_type, char, ::xsd::cxx::tree::schema_type::double_ > deltaT_traits;

  const deltaT_type&
  deltaT () const;

  deltaT_type&
  deltaT ();

  void
  deltaT (const deltaT_type& x);

  // inputType
  //
  typedef ::xml_schema::string inputType_type;
  typedef ::xsd::cxx::tree::traits< inputType_type, char > inputType_traits;

  const inputType_type&
  inputType () const;

  inputType_type&
  inputType ();

  void
  inputType (const inputType_type& x);

  void
  inputType (::std::unique_ptr< inputType_type > p);

  // particleContainerType
  //
  typedef ::xml_schema::string particleContainerType_type;
  typedef ::xsd::cxx::tree::traits< particleContainerType_type, char > particleContainerType_traits;

  const particleContainerType_type&
  particleContainerType () const;

  particleContainerType_type&
  particleContainerType ();

  void
  particleContainerType (const particleContainerType_type& x);

  void
  particleContainerType (::std::unique_ptr< particleContainerType_type > p);

  // r_cutoff
  //
  typedef ::xml_schema::double_ r_cutoff_type;
  typedef ::xsd::cxx::tree::traits< r_cutoff_type, char, ::xsd::cxx::tree::schema_type::double_ > r_cutoff_traits;

  const r_cutoff_type&
  r_cutoff () const;

  r_cutoff_type&
  r_cutoff ();

  void
  r_cutoff (const r_cutoff_type& x);

  // domainSizeX
  //
  typedef ::xml_schema::int_ domainSizeX_type;
  typedef ::xsd::cxx::tree::traits< domainSizeX_type, char > domainSizeX_traits;

  const domainSizeX_type&
  domainSizeX () const;

  domainSizeX_type&
  domainSizeX ();

  void
  domainSizeX (const domainSizeX_type& x);

  // domainSizeY
  //
  typedef ::xml_schema::int_ domainSizeY_type;
  typedef ::xsd::cxx::tree::traits< domainSizeY_type, char > domainSizeY_traits;

  const domainSizeY_type&
  domainSizeY () const;

  domainSizeY_type&
  domainSizeY ();

  void
  domainSizeY (const domainSizeY_type& x);

  // domainSizeZ
  //
  typedef ::xml_schema::int_ domainSizeZ_type;
  typedef ::xsd::cxx::tree::traits< domainSizeZ_type, char > domainSizeZ_traits;

  const domainSizeZ_type&
  domainSizeZ () const;

  domainSizeZ_type&
  domainSizeZ ();

  void
  domainSizeZ (const domainSizeZ_type& x);

  // boundary1Type
  //
  typedef ::xml_schema::int_ boundary1Type_type;
  typedef ::xsd::cxx::tree::traits< boundary1Type_type, char > boundary1Type_traits;

  const boundary1Type_type&
  boundary1Type () const;

  boundary1Type_type&
  boundary1Type ();

  void
  boundary1Type (const boundary1Type_type& x);

  // boundary2Type
  //
  typedef ::xml_schema::int_ boundary2Type_type;
  typedef ::xsd::cxx::tree::traits< boundary2Type_type, char > boundary2Type_traits;

  const boundary2Type_type&
  boundary2Type () const;

  boundary2Type_type&
  boundary2Type ();

  void
  boundary2Type (const boundary2Type_type& x);

  // boundary3Type
  //
  typedef ::xml_schema::int_ boundary3Type_type;
  typedef ::xsd::cxx::tree::traits< boundary3Type_type, char > boundary3Type_traits;

  const boundary3Type_type&
  boundary3Type () const;

  boundary3Type_type&
  boundary3Type ();

  void
  boundary3Type (const boundary3Type_type& x);

  // boundary4Type
  //
  typedef ::xml_schema::int_ boundary4Type_type;
  typedef ::xsd::cxx::tree::traits< boundary4Type_type, char > boundary4Type_traits;

  const boundary4Type_type&
  boundary4Type () const;

  boundary4Type_type&
  boundary4Type ();

  void
  boundary4Type (const boundary4Type_type& x);

  // boundary5Type
  //
  typedef ::xml_schema::int_ boundary5Type_type;
  typedef ::xsd::cxx::tree::traits< boundary5Type_type, char > boundary5Type_traits;

  const boundary5Type_type&
  boundary5Type () const;

  boundary5Type_type&
  boundary5Type ();

  void
  boundary5Type (const boundary5Type_type& x);

  // boundary6Type
  //
  typedef ::xml_schema::int_ boundary6Type_type;
  typedef ::xsd::cxx::tree::traits< boundary6Type_type, char > boundary6Type_traits;

  const boundary6Type_type&
  boundary6Type () const;

  boundary6Type_type&
  boundary6Type ();

  void
  boundary6Type (const boundary6Type_type& x);

  // temp_init
  //
  typedef ::xml_schema::double_ temp_init_type;
  typedef ::xsd::cxx::tree::traits< temp_init_type, char, ::xsd::cxx::tree::schema_type::double_ > temp_init_traits;

  const temp_init_type&
  temp_init () const;

  temp_init_type&
  temp_init ();

  void
  temp_init (const temp_init_type& x);

  // n_thermostat
  //
  typedef ::xml_schema::int_ n_thermostat_type;
  typedef ::xsd::cxx::tree::traits< n_thermostat_type, char > n_thermostat_traits;

  const n_thermostat_type&
  n_thermostat () const;

  n_thermostat_type&
  n_thermostat ();

  void
  n_thermostat (const n_thermostat_type& x);

  // temp_target
  //
  typedef ::xml_schema::double_ temp_target_type;
  typedef ::xsd::cxx::tree::traits< temp_target_type, char, ::xsd::cxx::tree::schema_type::double_ > temp_target_traits;

  const temp_target_type&
  temp_target () const;

  temp_target_type&
  temp_target ();

  void
  temp_target (const temp_target_type& x);

  // delta_temp
  //
  typedef ::xml_schema::double_ delta_temp_type;
  typedef ::xsd::cxx::tree::traits< delta_temp_type, char, ::xsd::cxx::tree::schema_type::double_ > delta_temp_traits;

  const delta_temp_type&
  delta_temp () const;

  delta_temp_type&
  delta_temp ();

  void
  delta_temp (const delta_temp_type& x);

  // g_grav
  //
  typedef ::xml_schema::double_ g_grav_type;
  typedef ::xsd::cxx::tree::traits< g_grav_type, char, ::xsd::cxx::tree::schema_type::double_ > g_grav_traits;

  const g_grav_type&
  g_grav () const;

  g_grav_type&
  g_grav ();

  void
  g_grav (const g_grav_type& x);

  // particles
  //
  typedef ::particles particles_type;
  typedef ::xsd::cxx::tree::sequence< particles_type > particles_sequence;
  typedef particles_sequence::iterator particles_iterator;
  typedef particles_sequence::const_iterator particles_const_iterator;
  typedef ::xsd::cxx::tree::traits< particles_type, char > particles_traits;

  const particles_sequence&
  particles () const;

  particles_sequence&
  particles ();

  void
  particles (const particles_sequence& s);

  // cuboids
  //
  typedef ::cuboids cuboids_type;
  typedef ::xsd::cxx::tree::sequence< cuboids_type > cuboids_sequence;
  typedef cuboids_sequence::iterator cuboids_iterator;
  typedef cuboids_sequence::const_iterator cuboids_const_iterator;
  typedef ::xsd::cxx::tree::traits< cuboids_type, char > cuboids_traits;

  const cuboids_sequence&
  cuboids () const;

  cuboids_sequence&
  cuboids ();

  void
  cuboids (const cuboids_sequence& s);

  // disk
  //
  typedef ::disk disk_type;
  typedef ::xsd::cxx::tree::sequence< disk_type > disk_sequence;
  typedef disk_sequence::iterator disk_iterator;
  typedef disk_sequence::const_iterator disk_const_iterator;
  typedef ::xsd::cxx::tree::traits< disk_type, char > disk_traits;

  const disk_sequence&
  disk () const;

  disk_sequence&
  disk ();

  void
  disk (const disk_sequence& s);

  // Constructors.
  //
  input (const tStart_type&,
         const tEnd_type&,
         const deltaT_type&,
         const inputType_type&,
         const particleContainerType_type&,
         const r_cutoff_type&,
         const domainSizeX_type&,
         const domainSizeY_type&,
         const domainSizeZ_type&,
         const boundary1Type_type&,
         const boundary2Type_type&,
         const boundary3Type_type&,
         const boundary4Type_type&,
         const boundary5Type_type&,
         const boundary6Type_type&,
         const temp_init_type&,
         const n_thermostat_type&,
         const temp_target_type&,
         const delta_temp_type&,
         const g_grav_type&);

  input (const ::xercesc::DOMElement& e,
         ::xml_schema::flags f = 0,
         ::xml_schema::container* c = 0);

  input (const input& x,
         ::xml_schema::flags f = 0,
         ::xml_schema::container* c = 0);

  virtual input*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  input&
  operator= (const input& x);

  virtual 
  ~input ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< tStart_type > tStart_;
  ::xsd::cxx::tree::one< tEnd_type > tEnd_;
  ::xsd::cxx::tree::one< deltaT_type > deltaT_;
  ::xsd::cxx::tree::one< inputType_type > inputType_;
  ::xsd::cxx::tree::one< particleContainerType_type > particleContainerType_;
  ::xsd::cxx::tree::one< r_cutoff_type > r_cutoff_;
  ::xsd::cxx::tree::one< domainSizeX_type > domainSizeX_;
  ::xsd::cxx::tree::one< domainSizeY_type > domainSizeY_;
  ::xsd::cxx::tree::one< domainSizeZ_type > domainSizeZ_;
  ::xsd::cxx::tree::one< boundary1Type_type > boundary1Type_;
  ::xsd::cxx::tree::one< boundary2Type_type > boundary2Type_;
  ::xsd::cxx::tree::one< boundary3Type_type > boundary3Type_;
  ::xsd::cxx::tree::one< boundary4Type_type > boundary4Type_;
  ::xsd::cxx::tree::one< boundary5Type_type > boundary5Type_;
  ::xsd::cxx::tree::one< boundary6Type_type > boundary6Type_;
  ::xsd::cxx::tree::one< temp_init_type > temp_init_;
  ::xsd::cxx::tree::one< n_thermostat_type > n_thermostat_;
  ::xsd::cxx::tree::one< temp_target_type > temp_target_;
  ::xsd::cxx::tree::one< delta_temp_type > delta_temp_;
  ::xsd::cxx::tree::one< g_grav_type > g_grav_;
  particles_sequence particles_;
  cuboids_sequence cuboids_;
  disk_sequence disk_;
};

class output: public ::xml_schema::type
{
  public:
  // baseName
  //
  typedef ::xml_schema::string baseName_type;
  typedef ::xsd::cxx::tree::traits< baseName_type, char > baseName_traits;

  const baseName_type&
  baseName () const;

  baseName_type&
  baseName ();

  void
  baseName (const baseName_type& x);

  void
  baseName (::std::unique_ptr< baseName_type > p);

  // writeFrequency
  //
  typedef ::xml_schema::int_ writeFrequency_type;
  typedef ::xsd::cxx::tree::traits< writeFrequency_type, char > writeFrequency_traits;

  const writeFrequency_type&
  writeFrequency () const;

  writeFrequency_type&
  writeFrequency ();

  void
  writeFrequency (const writeFrequency_type& x);

  // outputType
  //
  typedef ::xml_schema::string outputType_type;
  typedef ::xsd::cxx::tree::traits< outputType_type, char > outputType_traits;

  const outputType_type&
  outputType () const;

  outputType_type&
  outputType ();

  void
  outputType (const outputType_type& x);

  void
  outputType (::std::unique_ptr< outputType_type > p);

  // Constructors.
  //
  output (const baseName_type&,
          const writeFrequency_type&,
          const outputType_type&);

  output (const ::xercesc::DOMElement& e,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  output (const output& x,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  virtual output*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  output&
  operator= (const output& x);

  virtual 
  ~output ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< baseName_type > baseName_;
  ::xsd::cxx::tree::one< writeFrequency_type > writeFrequency_;
  ::xsd::cxx::tree::one< outputType_type > outputType_;
};

class config: public ::xml_schema::type
{
  public:
  // performanceMeasurement
  //
  typedef ::xml_schema::boolean performanceMeasurement_type;
  typedef ::xsd::cxx::tree::optional< performanceMeasurement_type > performanceMeasurement_optional;
  typedef ::xsd::cxx::tree::traits< performanceMeasurement_type, char > performanceMeasurement_traits;

  const performanceMeasurement_optional&
  performanceMeasurement () const;

  performanceMeasurement_optional&
  performanceMeasurement ();

  void
  performanceMeasurement (const performanceMeasurement_type& x);

  void
  performanceMeasurement (const performanceMeasurement_optional& x);

  // logLevel
  //
  typedef ::xml_schema::string logLevel_type;
  typedef ::xsd::cxx::tree::optional< logLevel_type > logLevel_optional;
  typedef ::xsd::cxx::tree::traits< logLevel_type, char > logLevel_traits;

  const logLevel_optional&
  logLevel () const;

  logLevel_optional&
  logLevel ();

  void
  logLevel (const logLevel_type& x);

  void
  logLevel (const logLevel_optional& x);

  void
  logLevel (::std::unique_ptr< logLevel_type > p);

  // Constructors.
  //
  config ();

  config (const ::xercesc::DOMElement& e,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  config (const config& x,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  virtual config*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  config&
  operator= (const config& x);

  virtual 
  ~config ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  performanceMeasurement_optional performanceMeasurement_;
  logLevel_optional logLevel_;
};

class particles: public ::xml_schema::type
{
  public:
  // x
  //
  typedef ::xml_schema::double_ x_type;
  typedef ::xsd::cxx::tree::traits< x_type, char, ::xsd::cxx::tree::schema_type::double_ > x_traits;

  const x_type&
  x () const;

  x_type&
  x ();

  void
  x (const x_type& x);

  // y
  //
  typedef ::xml_schema::double_ y_type;
  typedef ::xsd::cxx::tree::traits< y_type, char, ::xsd::cxx::tree::schema_type::double_ > y_traits;

  const y_type&
  y () const;

  y_type&
  y ();

  void
  y (const y_type& x);

  // z
  //
  typedef ::xml_schema::double_ z_type;
  typedef ::xsd::cxx::tree::traits< z_type, char, ::xsd::cxx::tree::schema_type::double_ > z_traits;

  const z_type&
  z () const;

  z_type&
  z ();

  void
  z (const z_type& x);

  // velocityX
  //
  typedef ::xml_schema::double_ velocityX_type;
  typedef ::xsd::cxx::tree::traits< velocityX_type, char, ::xsd::cxx::tree::schema_type::double_ > velocityX_traits;

  const velocityX_type&
  velocityX () const;

  velocityX_type&
  velocityX ();

  void
  velocityX (const velocityX_type& x);

  // velocityY
  //
  typedef ::xml_schema::double_ velocityY_type;
  typedef ::xsd::cxx::tree::traits< velocityY_type, char, ::xsd::cxx::tree::schema_type::double_ > velocityY_traits;

  const velocityY_type&
  velocityY () const;

  velocityY_type&
  velocityY ();

  void
  velocityY (const velocityY_type& x);

  // velocityZ
  //
  typedef ::xml_schema::double_ velocityZ_type;
  typedef ::xsd::cxx::tree::traits< velocityZ_type, char, ::xsd::cxx::tree::schema_type::double_ > velocityZ_traits;

  const velocityZ_type&
  velocityZ () const;

  velocityZ_type&
  velocityZ ();

  void
  velocityZ (const velocityZ_type& x);

  // mass
  //
  typedef ::xml_schema::double_ mass_type;
  typedef ::xsd::cxx::tree::traits< mass_type, char, ::xsd::cxx::tree::schema_type::double_ > mass_traits;

  const mass_type&
  mass () const;

  mass_type&
  mass ();

  void
  mass (const mass_type& x);

  // Constructors.
  //
  particles (const x_type&,
             const y_type&,
             const z_type&,
             const velocityX_type&,
             const velocityY_type&,
             const velocityZ_type&,
             const mass_type&);

  particles (const ::xercesc::DOMElement& e,
             ::xml_schema::flags f = 0,
             ::xml_schema::container* c = 0);

  particles (const particles& x,
             ::xml_schema::flags f = 0,
             ::xml_schema::container* c = 0);

  virtual particles*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  particles&
  operator= (const particles& x);

  virtual 
  ~particles ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< x_type > x_;
  ::xsd::cxx::tree::one< y_type > y_;
  ::xsd::cxx::tree::one< z_type > z_;
  ::xsd::cxx::tree::one< velocityX_type > velocityX_;
  ::xsd::cxx::tree::one< velocityY_type > velocityY_;
  ::xsd::cxx::tree::one< velocityZ_type > velocityZ_;
  ::xsd::cxx::tree::one< mass_type > mass_;
};

class cuboids: public ::xml_schema::type
{
  public:
  // n1
  //
  typedef ::xml_schema::double_ n1_type;
  typedef ::xsd::cxx::tree::traits< n1_type, char, ::xsd::cxx::tree::schema_type::double_ > n1_traits;

  const n1_type&
  n1 () const;

  n1_type&
  n1 ();

  void
  n1 (const n1_type& x);

  // n2
  //
  typedef ::xml_schema::double_ n2_type;
  typedef ::xsd::cxx::tree::traits< n2_type, char, ::xsd::cxx::tree::schema_type::double_ > n2_traits;

  const n2_type&
  n2 () const;

  n2_type&
  n2 ();

  void
  n2 (const n2_type& x);

  // n3
  //
  typedef ::xml_schema::double_ n3_type;
  typedef ::xsd::cxx::tree::traits< n3_type, char, ::xsd::cxx::tree::schema_type::double_ > n3_traits;

  const n3_type&
  n3 () const;

  n3_type&
  n3 ();

  void
  n3 (const n3_type& x);

  // distance
  //
  typedef ::xml_schema::double_ distance_type;
  typedef ::xsd::cxx::tree::traits< distance_type, char, ::xsd::cxx::tree::schema_type::double_ > distance_traits;

  const distance_type&
  distance () const;

  distance_type&
  distance ();

  void
  distance (const distance_type& x);

  // meanVelocity
  //
  typedef ::xml_schema::double_ meanVelocity_type;
  typedef ::xsd::cxx::tree::traits< meanVelocity_type, char, ::xsd::cxx::tree::schema_type::double_ > meanVelocity_traits;

  const meanVelocity_type&
  meanVelocity () const;

  meanVelocity_type&
  meanVelocity ();

  void
  meanVelocity (const meanVelocity_type& x);

  // dimension
  //
  typedef ::xml_schema::double_ dimension_type;
  typedef ::xsd::cxx::tree::traits< dimension_type, char, ::xsd::cxx::tree::schema_type::double_ > dimension_traits;

  const dimension_type&
  dimension () const;

  dimension_type&
  dimension ();

  void
  dimension (const dimension_type& x);

  // epsilon
  //
  typedef ::xml_schema::double_ epsilon_type;
  typedef ::xsd::cxx::tree::traits< epsilon_type, char, ::xsd::cxx::tree::schema_type::double_ > epsilon_traits;

  const epsilon_type&
  epsilon () const;

  epsilon_type&
  epsilon ();

  void
  epsilon (const epsilon_type& x);

  // sigma
  //
  typedef ::xml_schema::double_ sigma_type;
  typedef ::xsd::cxx::tree::traits< sigma_type, char, ::xsd::cxx::tree::schema_type::double_ > sigma_traits;

  const sigma_type&
  sigma () const;

  sigma_type&
  sigma ();

  void
  sigma (const sigma_type& x);

  // Constructors.
  //
  cuboids (const n1_type&,
           const n2_type&,
           const n3_type&,
           const distance_type&,
           const meanVelocity_type&,
           const dimension_type&,
           const epsilon_type&,
           const sigma_type&);

  cuboids (const ::xercesc::DOMElement& e,
           ::xml_schema::flags f = 0,
           ::xml_schema::container* c = 0);

  cuboids (const cuboids& x,
           ::xml_schema::flags f = 0,
           ::xml_schema::container* c = 0);

  virtual cuboids*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  cuboids&
  operator= (const cuboids& x);

  virtual 
  ~cuboids ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< n1_type > n1_;
  ::xsd::cxx::tree::one< n2_type > n2_;
  ::xsd::cxx::tree::one< n3_type > n3_;
  ::xsd::cxx::tree::one< distance_type > distance_;
  ::xsd::cxx::tree::one< meanVelocity_type > meanVelocity_;
  ::xsd::cxx::tree::one< dimension_type > dimension_;
  ::xsd::cxx::tree::one< epsilon_type > epsilon_;
  ::xsd::cxx::tree::one< sigma_type > sigma_;
};

class disk: public ::xml_schema::type
{
  public:
  // radius
  //
  typedef ::xml_schema::double_ radius_type;
  typedef ::xsd::cxx::tree::traits< radius_type, char, ::xsd::cxx::tree::schema_type::double_ > radius_traits;

  const radius_type&
  radius () const;

  radius_type&
  radius ();

  void
  radius (const radius_type& x);

  // distance
  //
  typedef ::xml_schema::double_ distance_type;
  typedef ::xsd::cxx::tree::traits< distance_type, char, ::xsd::cxx::tree::schema_type::double_ > distance_traits;

  const distance_type&
  distance () const;

  distance_type&
  distance ();

  void
  distance (const distance_type& x);

  // meanVelocity
  //
  typedef ::xml_schema::double_ meanVelocity_type;
  typedef ::xsd::cxx::tree::traits< meanVelocity_type, char, ::xsd::cxx::tree::schema_type::double_ > meanVelocity_traits;

  const meanVelocity_type&
  meanVelocity () const;

  meanVelocity_type&
  meanVelocity ();

  void
  meanVelocity (const meanVelocity_type& x);

  // dimension
  //
  typedef ::xml_schema::double_ dimension_type;
  typedef ::xsd::cxx::tree::traits< dimension_type, char, ::xsd::cxx::tree::schema_type::double_ > dimension_traits;

  const dimension_type&
  dimension () const;

  dimension_type&
  dimension ();

  void
  dimension (const dimension_type& x);

  // epsilon
  //
  typedef ::xml_schema::double_ epsilon_type;
  typedef ::xsd::cxx::tree::traits< epsilon_type, char, ::xsd::cxx::tree::schema_type::double_ > epsilon_traits;

  const epsilon_type&
  epsilon () const;

  epsilon_type&
  epsilon ();

  void
  epsilon (const epsilon_type& x);

  // sigma
  //
  typedef ::xml_schema::double_ sigma_type;
  typedef ::xsd::cxx::tree::traits< sigma_type, char, ::xsd::cxx::tree::schema_type::double_ > sigma_traits;

  const sigma_type&
  sigma () const;

  sigma_type&
  sigma ();

  void
  sigma (const sigma_type& x);

  // Constructors.
  //
  disk (const radius_type&,
        const distance_type&,
        const meanVelocity_type&,
        const dimension_type&,
        const epsilon_type&,
        const sigma_type&);

  disk (const ::xercesc::DOMElement& e,
        ::xml_schema::flags f = 0,
        ::xml_schema::container* c = 0);

  disk (const disk& x,
        ::xml_schema::flags f = 0,
        ::xml_schema::container* c = 0);

  virtual disk*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  disk&
  operator= (const disk& x);

  virtual 
  ~disk ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< radius_type > radius_;
  ::xsd::cxx::tree::one< distance_type > distance_;
  ::xsd::cxx::tree::one< meanVelocity_type > meanVelocity_;
  ::xsd::cxx::tree::one< dimension_type > dimension_;
  ::xsd::cxx::tree::one< epsilon_type > epsilon_;
  ::xsd::cxx::tree::one< sigma_type > sigma_;
};

#include <iosfwd>

#include <xercesc/sax/InputSource.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMErrorHandler.hpp>

// Parse a URI or a local file.
//

::std::unique_ptr< ::simulation >
simulation_ (const ::std::string& uri,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::simulation >
simulation_ (const ::std::string& uri,
             ::xml_schema::error_handler& eh,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::simulation >
simulation_ (const ::std::string& uri,
             ::xercesc::DOMErrorHandler& eh,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse std::istream.
//

::std::unique_ptr< ::simulation >
simulation_ (::std::istream& is,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::simulation >
simulation_ (::std::istream& is,
             ::xml_schema::error_handler& eh,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::simulation >
simulation_ (::std::istream& is,
             ::xercesc::DOMErrorHandler& eh,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::simulation >
simulation_ (::std::istream& is,
             const ::std::string& id,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::simulation >
simulation_ (::std::istream& is,
             const ::std::string& id,
             ::xml_schema::error_handler& eh,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::simulation >
simulation_ (::std::istream& is,
             const ::std::string& id,
             ::xercesc::DOMErrorHandler& eh,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse xercesc::InputSource.
//

::std::unique_ptr< ::simulation >
simulation_ (::xercesc::InputSource& is,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::simulation >
simulation_ (::xercesc::InputSource& is,
             ::xml_schema::error_handler& eh,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::simulation >
simulation_ (::xercesc::InputSource& is,
             ::xercesc::DOMErrorHandler& eh,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse xercesc::DOMDocument.
//

::std::unique_ptr< ::simulation >
simulation_ (const ::xercesc::DOMDocument& d,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::simulation >
simulation_ (::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d,
             ::xml_schema::flags f = 0,
             const ::xml_schema::properties& p = ::xml_schema::properties ());

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

#endif // SIMULATION_HXX
