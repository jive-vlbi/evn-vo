<resource schema="evn">
  <meta name="title">EVN archive</meta>
  <meta name="coverage">
    <meta name="profile">AllSky J2000</meta>
    <meta name="waveband">Radio</meta>
  </meta>

  <table id="main" onDisk="True" namePath="//obscore#ObsCore">
    <meta name="description">
      Visibility data from the EVN archive.
    </meta>
    <meta name="_associatedDatalinkService">
      <meta name="serviceId">dl</meta>
      <meta name="idColumn">obs_publisher_did</meta>
    </meta>
    <mixin>//products#table</mixin>
    <mixin
	calibLevel="1"
	collectionName="'EVN'"
	instrumentName="'EVN'"
	oUCD="'stat.uncalib'"
	productType="'visibility'"
	sXel1="-1"
	sXel2="-1"
	accessURL="access_url"
	coverage="s_region"
	dec="s_dec"
	did="obs_publisher_did"
	emMax="em_max"
	emMin="em_min"
	emResPower="em_res_power"
	emXel="em_xel"
	expTime="t_exptime"
	fov="s_fov"
	obsId="obs_id"
	polXel="pol_xel"
	ra="s_ra"
	size="access_estsize"
	sResolution="s_resolution"
	targetName="target_name"
	tMax="t_max"
	tMin="t_min"
	tResolution="t_resolution"
	tXel="t_xel"
	preview="NULL"
	>//obscore#publish</mixin>

    <column name="_nparts" type="smallint" required="True"
	    ucd="meta.number" verbLevel="30"/>

    <column original="accref" displayHint=""/>

    <!-- Need to use the obscore-extracolumns in useronfig.rd to achieve the
	 same thing in the ivoa.obscore table. -->
    <column original="access_url" displayHint="type=url">
	    <property name="anchorText">datalink</property>
    </column>

    <LOOP listItems="access_url s_region s_dec obs_publisher_did em_max em_min em_res_power em_xel t_exptime s_fov obs_id pol_xel s_ra access_estsize s_resolution target_name t_max t_min t_resolution t_xel">
      <events>
	<column original="\item"/>
      </events>
    </LOOP>

  </table>
  
  <data id="import">
    <sources>data/evn.csv</sources>
    <csvGrammar>
      <rowfilter name="makeAccref">
	<setup>
	  <code>
	    import urllib
	  </code>
	</setup>
	<code>
	  @accref = getAccrefFromStandardPubDID(@obs_publisher_did)
	  @pubDID = @obs_publisher_did
	  yield row
	</code>
      </rowfilter>

      <rowfilter procDef="//products#define">
	<bind key="table">"evn.main"</bind>
	<bind key="accref">@accref</bind>
	<bind key="mime">"application/x-votable+xml;content=datalink"</bind>
      </rowfilter>
    </csvGrammar>
    <make table="main">
      <rowmaker idmaps="*">
	<var key="s_region">pgsphere.SCircle(pgsphere.SPoint.fromDegrees(float(@s_ra), float(@s_dec)), float(@s_fov)).asPoly()</var>
	<var key="access_url">makeAbsoluteURL(
		"\rdId/dl/dlmeta?ID="+urllib.parse.quote(@pubDID))</var>
      </rowmaker>
    </make>
  </data>

  <service id="dl" allowed="dlmeta">
    <meta name="title">EVN Visibility Datalink Service</meta>
    <datalinkCore>
      <metaMaker>
	<code>
	  with base.getTableConn() as conn:
	    query = "SELECT * from %s WHERE obs_publisher_did='%s'" % (descriptor.sourceTable, descriptor.pubDID)
	    for row in conn.queryToDicts(query):
	      nparts=row['_nparts']
	      break
	  accref = getAccrefFromStandardPubDID(descriptor.pubDID)
	  pos = accref.find('_', accref.find('_') + 1)
	  exp_id = accref[0:pos]
	  product_id = accref[pos:]
	  pos = accref.find('_')
	  exp = accref[0:pos]
	  baseurl = "http://archive.jive.nl/exp/" + exp_id + "/fits/" + exp.lower() + product_id + ".IDI"
	  descriptor.suppressAutoLinks = True
	  for part in range(nparts):
	    url = baseurl + "%d" % (part + 1)
	    yield LinkDef(descriptor.pubDID, url, semantics="#this", description="FITS-IDI", contentType="application/x-fits-idi")
	</code>
      </metaMaker>
    </datalinkCore>
  </service>

</resource>
