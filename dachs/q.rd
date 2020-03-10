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
    <mixin>//products#table</mixin>
    <mixin
	calibLevel="1"
	collectionName="'EVN'"
	instrumentName="'EVN'"
	mime="'application/x-fits-idi'"
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
	>//obscore#publish</mixin>
    <LOOP listItems="access_url s_region s_dec obs_publisher_did em_max em_min em_res_power em_xel t_exptime s_fov obs_id pol_xel s_ra access_estsize s_resolution target_name t_max t_min t_resolution t_xel">
      <events>
	<column original="\item"/>
      </events>
    </LOOP>
  </table>
  
  <data id="import">
    <sources>data/evn.csv</sources>
    <csvGrammar>
      <rowfilter procDef="//products#define">
	<bind key="table">"evn.main"</bind>
      </rowfilter>
    </csvGrammar>
    <make table="main">
      <rowmaker idmaps="*">
	<var key="s_region">pgsphere.SCircle(pgsphere.SPoint.fromDegrees(float(@s_ra), float(@s_dec)), float(@s_fov)).asPoly()</var>
      </rowmaker>
    </make>
  </data>

  <service id="vis">
    <meta name="title">EVN Visibility Data</meta>
    <dbCore queriedTable="main">
    </dbCore>
  </service>

</resource>
