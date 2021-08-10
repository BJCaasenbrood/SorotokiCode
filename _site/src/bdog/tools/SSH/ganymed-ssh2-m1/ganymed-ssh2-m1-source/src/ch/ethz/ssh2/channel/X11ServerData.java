
package ch.ethz.ssh2.channel;

/**
 * X11ServerData. Data regarding an x11 forwarding target.
 *
 * @author Christian Plattner
 * @version 2.50, 03/15/10
 * 
 */
public class X11ServerData
{
	public String hostname;
	public int port;
	public byte[] x11_magic_cookie; /* not the remote (fake) one, the local (real) one */
}
