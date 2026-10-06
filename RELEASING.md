# Releasing jOSQP

Releases go to Maven Central through the Sonatype Central Portal. Everything is
driven by the `release` profile in `pom.xml`; the only things that live outside
the repository are the GPG key and the Portal credentials. **Never commit those,
and never paste them into a ticket or a chat.**

## One-time setup

1. **GPG key** for signing (Central rejects unsigned artifacts):

   ```
   gpg --full-generate-key          # RSA 4096, a real name and e-mail, no expiry is fine
   gpg --list-secret-keys --keyid-format=long
   gpg --keyserver keyserver.ubuntu.com --send-keys <KEY-ID>
   ```

2. **Central Portal account** at <https://central.sonatype.com>. The namespace
   `com.quantego` is already verified - `com.quantego:josqp:0.6.5` was published
   from it. A second namespace needs a DNS TXT record or a GitHub verification.

3. **User token**: Portal -> Account -> *Generate User Token*. It yields a
   username and a password; that pair is the `central` server credential.

4. **`~/.m2/settings.xml`** with the token (file permissions `600`):

   ```xml
   <settings>
     <servers>
       <server>
         <id>central</id>
         <username>TOKEN-USERNAME</username>
         <password>TOKEN-PASSWORD</password>
       </server>
     </servers>
   </settings>
   ```

5. **Maven and a JDK** on the PATH. The build targets Java 8 release level, so
   any JDK 9+ works; JDK 26 was used for the 1.0.0 verification.

## Release

```bash
mvn -Prelease -DskipTests deploy
```

The profile attaches the sources and javadoc jars, signs all four artifacts,
and uploads them as a *staged* deployment (`autoPublish=false`). The deployment
then has to be published explicitly in the Portal under *Deployments*; if
anything looks wrong, drop it there instead. Nothing is public before that
step, and a dropped deployment can be repeated with the same version.

The passphrase comes from `MAVEN_GPG_PASSPHRASE`, `MAVEN_GPG_PASSPHRASE_FILE`
or a running `gpg-agent`.

Do **not** add `-Pvector` here: the Vector API backend needs JDK 26 and the
`jdk.incubator.vector` module and must not end up in the plain Java 8 artifact.

## Afterwards

```bash
git tag -a v1.0.0 -m "jOSQP 1.0.0"
git push origin master --tags
```

Then create the GitHub release from the tag and copy the matching
[CHANGELOG](CHANGELOG.md) section as the release notes, and open a new
`## [Unreleased]` section for the next cycle.

## Verifying a release locally

```bash
mvn -Prelease -DskipTests package   # main, sources and javadoc jars in target/
mvn test                            # 90 tests
mvn -Pvector test                   # 93 tests, JDK 26
```
