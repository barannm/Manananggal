# Manananggal

Manananggal is a JAVA application developed to facilitate the identification of alternative splicing events (ASEs). The web application allows the visual inspection of hundreds of RNAseq samples to identify ASEs. The Manananggal assists you with this by generating a list of potential ASEs that you can review and share with others.
The source code necessary to compile the program and pre-compiled binary files are available in this repository.

Compiled program<br>
<a href=https://github.com/barannm/Manananggal/tree/master/war_files>Web Application</a>(.war file)<br>
<a href=https://github.com/barannm/Manananggal/tree/master/jar_files>Console Application</a>(.jar file)

If you are unfamilar with the use of .war and .jar files, please review how to deploy .WAR files in the user guide of your server installation (e.g. for a <a href=https://tomcat.apache.org/tomcat-6.0-doc/deployer-howto.html>Apache Tomcat Server</a> and the installation section in the <a href=https://services.bio.ifi.lmu.de/manananggal/Manual/index.html>Manananggal User Manual</a>, which is also available for download <a href=https://github.com/barannm/Manananggal/tree/master/WebContent/Manual>here</a>.

For the compilation of the source code several external libraries are required. These libraries include the <a href=https://commons.apache.org/>commons-lang3 and commons-math3</a> libraries (<a href=http://commons.apache.org/proper/commons-daemon/license.html>license</a>), the <a href=http://www.slf4j.org/>Simple Logging Facade for Java (SLF4J)</a> (<a href=http://www.slf4j.org/license.html>license</a>) library, <a href=https://www.broadinstitute.org/igv/igvtools>IGVTools</a> (<a href=https://opensource.org/licenses/lgpl-2.1.php>LGPL license</a>) and the <a href=https://www.zkoss.org/>ZK Community Edition</a> (<a href=https://opensource.org/licenses/lgpl-2.1.php>LGPL license</a>). and can thus be downloaded free of charge. Versions of these libraries used for the tool are included in the WebContent/WEB-INF/lib folder.
